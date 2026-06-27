import os
import zipfile
import xml.etree.ElementTree as ET
import sys
from collections import defaultdict

ASY_HEADER = """import geometry;
unitsize(1cm);
defaultpen(fontsize(10pt));

// --- 0. STYLE & COLOR DEFINITIONS ---
// DYNAMIC_STYLES_HERE

// --- GEO2D GGB-TO-ASY HELPER FUNCTIONS ---
// DYNAMIC_HELPERS_HERE
// -----------------------------------------
"""

ASY_HELPERS = {
    'ray': "line ray(point A, point B) { return line(A, false, B, true); }",
    'circle_point': "circle circle_point(point C, point P) { return circle(C, abs(C - P)); }",
    
    'triangle_center': """point triangle_center(point A, point B, point C, int center_type) {
    if (center_type == 1001) return excenter(A, B, C); // EXCENTER_A
    if (center_type == 1002) return excenter(B, C, A); // EXCENTER_B
    if (center_type == 1003) return excenter(C, A, B); // EXCENTER_C
    
    real weight(real a, real b, real c) {
        if (center_type == 1) return a; // INCENTER
        if (center_type == 2) return 1.0; // CENTROID
        if (center_type == 3) return a^2 * (b^2 + c^2 - a^2); // CIRCUMCENTER
        if (center_type == 4) return (a^2 + b^2 - c^2) * (a^2 + c^2 - b^2); // ORTHOCENTER
        if (center_type == 5) return a^2 * (b^2 + c^2) - (b^2 - c^2)^2; // NINE_POINT
        if (center_type == 6) return a^2; // SYMMEDIAN
        if (center_type == 7) return (b + c - a != 0) ? 1.0 / (b + c - a) : 0.0; // GERGONNE
        if (center_type == 8) return b + c - a; // NAGEL
        if (center_type == 9) return a * (b + c - a); // MITTENPUNKT
        if (center_type == 10) return b + c; // SPIEKER
        if (center_type == 11) return (b + c - a) * (b - c)^2; // FEUERBACH
        if (center_type == 69) return b^2 + c^2 - a^2; // X69
        if (center_type == 76) return (a != 0) ? 1.0 / a^2 : 0.0; // X76
        return 0.0;
    }

    real a = length(B - C);
    real b = length(C - A);
    real c = length(A - B);
    real wA = weight(a, b, c);
    real wB = weight(b, c, a);
    real wC = weight(c, a, b);
    real total = wA + wB + wC;
    if (abs(total) < 1e-12) return (nan, nan);
    return (wA*A + wB*B + wC*C) / total;
}""",

    'mixtilinear_incircle': """circle mixtilinear_incircle(point A, point B, point C) {
    real a = length(B - C); real b = length(C - A); real c = length(A - B);
    real s = (a + b + c) / 2.0;
    real area = sqrt(s * (s - a) * (s - b) * (s - c));
    real r = area / s;
    real cos_A = (b^2 + c^2 - a^2) / (2.0 * b * c);
    real radius = 2.0 * r / (1.0 + cos_A);
    point I = triangle_center(A, B, C, 1);
    point center = A + (radius / r) * (I - A);
    return circle(center, radius);
}""",

    'exmixtilinear_incircle': """circle exmixtilinear_incircle(point A, point B, point C) {
    real a = length(B - C); real b = length(C - A); real c = length(A - B);
    real s = (a + b + c) / 2.0;
    real area = sqrt(s * (s - a) * (s - b) * (s - c));
    real r_a = area / (s - a);
    real cos_A = (b^2 + c^2 - a^2) / (2.0 * b * c);
    real radius = 2.0 * r_a / (1.0 + cos_A);
    point Ia = triangle_center(A, B, C, 1001);
    point center = A + (radius / r_a) * (Ia - A);
    return circle(center, radius);
}""",

    'intersect_ggb': """// Helper for non-linear intersections with stable sorting by coordinates
point intersect_ggb(point[] pts, int idx) {
    if (pts.length == 0) return (nan, nan);
    // Stable bubble sort by X coordinate, then Y coordinate to match GGB standard indices
    for (int i = 0; i < pts.length - 1; ++i) {
        for (int j = i + 1; j < pts.length; ++j) {
            if (pts[i].x > pts[j].x || (abs(pts[i].x - pts[j].x) < 1e-9 && pts[i].y > pts[j].y)) {
                point temp = pts[i];
                pts[i] = pts[j];
                pts[j] = temp;
            }
        }
    }
    if (idx >= 1 && idx <= pts.length) return pts[idx-1];
    return pts[0];
}
point intersect_ggb(line l, circle c, int idx) { return intersect_ggb(intersectionpoints(l, c), idx); }
point intersect_ggb(circle c, line l, int idx) { return intersect_ggb(intersectionpoints(l, c), idx); }
point intersect_ggb(circle c1, circle c2, int idx) { return intersect_ggb(intersectionpoints(c1, c2), idx); }
point intersect_ggb(path p1, path p2, int idx) { return intersect_ggb(intersectionpoints(p1, p2), idx); }""",

    'relpoint_ggb_line': """// Project a static point onto a line
point relpoint_ggb_line(line l, point P) {
    return projection(l) * P;
}""",

    'relpoint_ggb_circle': """// Project a static point onto a circle analytically
point relpoint_ggb_circle(circle c, point P) {
    return c.C + c.r * unit(P - c.C);
}""",

    'relpoint_ggb_path': """// Project a static point onto a path using fast Ternary Search
point relpoint_ggb_path(path p, point P) {
    real L = 0;
    real R = length(p);
    for (int i=0; i<25; ++i) {
        real m1 = L + (R - L) / 3;
        real m2 = R - (R - L) / 3;
        real d1 = length(relpoint(p, m1) - P);
        real d2 = length(relpoint(p, m2) - P);
        if (d1 < d2)
            R = m2;
        else
            L = m1;
    }
    return relpoint(p, (L + R) / 2.0);
}""",

    'auto_align': """// Auto-align labels by pushing them radially outward from the centroid
pair auto_align(point P, point[] all_pts) {
    point centroid = (0,0);
    for (point pt : all_pts) centroid += pt;
    if (all_pts.length > 0) centroid /= all_pts.length;
    pair d = P - centroid;
    if (abs(d) < 1e-5) return NE;
    return unit(d);
}""",

    'polar_ggb': """// Polar line of a point with respect to a conic
line polar_ggb(point P, conic co) {
    bqe eq = equation(co);
    real x0 = P.x, y0 = P.y;
    real A = eq.a[0] * x0 + eq.a[1]/2 * y0 + eq.a[3]/2;
    real B = eq.a[1]/2 * x0 + eq.a[2] * y0 + eq.a[4]/2;
    real C = eq.a[3]/2 * x0 + eq.a[4]/2 * y0 + eq.a[5];
    return line(A, B, C);
}""",

    'arc_ggb': """// Circular arc through 3 points A, B, C
arc arc_ggb(point A, point B, point C) {
    circle c = circle(A, B, C);
    real cross_prod = (B.x - A.x)*(C.y - A.y) - (B.y - A.y)*(C.x - A.x);
    return arc(c, A, C, cross_prod > 0 ? CCW : CW);
}""",

    'tangents_ggb': """// Common tangents to two conics using dual conic intersection
real [][] tomatrix_ggb(conic c) {
    real []a = equation(c).a;
    real [][] m={{a[0], a[1]/2, a[3]/2},{a[1]/2, a[2], a[4]/2},{a[3]/2, a[4]/2, a[5]}};
    return m;
}

real [][] adjmatrix_ggb(real [][] m) {
    real [][] adj= {{m[1][1]*m[2][2] - m[1][2]*m[2][1], m[0][2]*m[2][1]- m[0][1]*m[2][2], m[0][1]*m[1][2]-m[0][2]*m[1][1]}, {m[1][2]*m[2][0] - m[1][0]*m[2][2], m[0][0]*m[2][2]- m[0][2]*m[2][0], m[0][2]*m[1][0]-m[0][0]*m[1][2]}, {m[1][0]*m[2][1] - m[1][1]*m[2][0], m[0][1]*m[2][0]- m[0][0]*m[2][1], m[0][0]*m[1][1]-m[0][1]*m[1][0]}};
    return adj;
}

real [] tocoefs_ggb(real [][] z) {
    real [] a = {z[0][0],2*z[0][1],z[1][1],2*z[0][2],2*z[1][2],z[2][2]};
    return a;
}

line [] tangents_ggb(conic cm, conic cn) {
    line [] tmp;
    real [][] m = tomatrix_ggb(cm), n = tomatrix_ggb(cn);
    real [] mc = tocoefs_ggb(adjmatrix_ggb(m)), nc=tocoefs_ggb(adjmatrix_ggb(n));
    bqe bm = bqe(mc[0],mc[1],mc[2],mc[3],mc[4],mc[5]);
    bqe bn = bqe(nc[0],nc[1],nc[2],nc[3],nc[4],nc[5]);
    point []p = intersectionpoints(conic(bm),conic(bn));
    for (int i=0; i<p.length; ++i) tmp.push(line(p[i].x, p[i].y, 1));
    return tmp;
}

// Common tangents to two circles using homothetic centers
line[] tangents_circle_circle(circle c1, circle c2) {
    line[] tmp;
    if (abs(c1.r - c2.r) > 1e-7) {
        point H_ext = (c1.r * c2.C - c2.r * c1.C) / (c1.r - c2.r);
        line[] t_ext = tangents(c1, H_ext);
        for (int i=0; i<t_ext.length; ++i) tmp.push(t_ext[i]);
    } else {
        vector v = unit(c2.C - c1.C);
        vector n = rotate(90)*v;
        tmp.push(line(c1.C + c1.r * n, c2.C + c2.r * n));
        tmp.push(line(c1.C - c1.r * n, c2.C - c2.r * n));
    }
    
    point H_int = (c1.r * c2.C + c2.r * c1.C) / (c1.r + c2.r);
    line[] t_int = tangents(c1, H_int);
    for (int i=0; i<t_int.length; ++i) tmp.push(t_int[i]);
    
    return tmp;
}""",

    # --- CUSTOM TOOLS HELPERS ---
    'harmonic_conjugate': """point harmonic_conjugate(point A, point B, point C) {
    real r = abs(C - A) / abs(C - B);
    real t = r / (r + 1.0);
    point D_int = (1.0 - t) * A + t * B;
    if (abs(C - D_int) < 1e-5) {
        real t_ext = (r != 1.0) ? r / (r - 1.0) : 1e5;
        return (1.0 - t_ext) * A + t_ext * B;
    } else {
        return D_int;
    }
}""",

    'circle_center_to_circle': """circle circle_center_to_circle(point O, circle c) {
    real d = distance(O, c.C);
    real r = abs(d - c.r);
    return circle(O, r);
}""",

    'locus_pl': """conic locus_pl(point P, line l) { return conic(P, l); } // Placeholder for custom conic locus""",
    'locus_pc': """conic locus_pc(point P, circle c) { return conic(P, c.C); } // Placeholder for custom conic locus""",
    'locus_lc': """conic locus_lc(line l, circle c) { return conic(c.C, l); } // Placeholder for custom conic locus""",
    'locus_cc': """conic locus_cc(circle c1, circle c2) { return conic(c1.C, c2.C); } // Placeholder for custom conic locus"""
}

def map_transform(ins, out, types, trans_type):
    target = ins[0]
    target_type = types.get(target, 'point')
    asy_type = 'line' if target_type == 'ray' else target_type
    transform = ""

    if trans_type == 'mirror':
        mirror_obj = ins[1]
        mirror_type = types.get(mirror_obj, 'point')
        
        if mirror_type == 'circle':
            if target_type in ['circle', 'line', 'segment', 'ray', 'conic']:
                return f"circle {out} = inversion({mirror_obj}) * {target};", "circle"
            else:
                return f"{asy_type} {out} = inversion({mirror_obj}) * {target};", target_type
        
        if mirror_type in ['line', 'segment', 'ray']:
            transform = f"reflect({mirror_obj})"
        else:
            transform = f"scale(-1, {mirror_obj})"
    
    elif trans_type == 'translate':
        transform = f"shift({ins[1]})"
        
    elif trans_type == 'rotate':
        raw_angle = str(ins[1])
        if '°' in raw_angle:
            angle = raw_angle.replace('°', '')
        else:
            angle = f"({raw_angle}) * 180 / pi"
            
        if len(ins) == 3:
            transform = f"rotate({angle}, {ins[2]})"
        else:
            transform = f"rotate({angle})"
            
    elif trans_type == 'dilate':
        if len(ins) == 3:
            transform = f"scale({ins[1]}, {ins[2]})"
        else:
            transform = f"scale({ins[1]})"

    if target_type == 'circle':
        if trans_type == 'dilate':
            return f"circle {out} = circle({transform} * {target}.C, abs({ins[1]}) * {target}.r);", "circle"
        else:
            return f"circle {out} = circle({transform} * {target}.C, {target}.r);", "circle"
    else:
        return f"{asy_type} {out} = {transform} * {target};", target_type


# GGB-TO-ASY COMMAND MAPPER GROUPS

# 1. POINT COMMANDS
POINT_COMMANDS = {
    'point': lambda ins, out, types: (f"point {out} = point({ins[0]});", "point"),
    'center': lambda ins, out, types: (f"point {out} = {ins[0]}.C;", "point"),
    'midpoint': lambda ins, out, types: (
        f"point {out} = midpoint(segment({ins[0]}, {ins[1]}));" if len(ins) == 2 else f"point {out} = {ins[0]}.C;", "point"
    ),
    'intersect': lambda ins, out, types: (
        f"point {out} = intersectionpoint({ins[0]}, {ins[1]});" if types.get(ins[0]) in ['line', 'segment', 'ray'] and types.get(ins[1]) in ['line', 'segment', 'ray'] else
        f"point {out} = intersect_ggb({ins[0]}, {ins[1]}, {int(ins[2]) if len(ins) >= 3 else 1});", "point"
    ),
    'intersection': lambda ins, out, types: (
        f"point {out} = intersectionpoint({ins[0]}, {ins[1]});" if types.get(ins[0]) in ['line', 'segment', 'ray'] and types.get(ins[1]) in ['line', 'segment', 'ray'] else
        f"point {out} = intersect_ggb({ins[0]}, {ins[1]}, {int(ins[2]) if len(ins) >= 3 else 1});", "point"
    ),
    'trianglecenter': lambda ins, out, types: (f"point {out} = triangle_center({ins[0]}, {ins[1]}, {ins[2]}, {ins[3]});", "point"),
    'closestpoint': lambda ins, out, types: (
        f"point {out} = relpoint_ggb_{'circle' if types.get(ins[0]) == 'circle' else ('line' if types.get(ins[0]) in ['line', 'segment', 'ray', 'axis'] else 'path')}({ins[0]}, {ins[1]});", "point"
    ),
}

# 2. LINE COMMANDS
LINE_COMMANDS = {
    'segment': lambda ins, out, types: (f"segment {out} = segment({ins[0]}, {ins[1]});", "segment"),
    'line': lambda ins, out, types: (
        f"line {out} = parallel({ins[0]}, {ins[1]});" if types.get(ins[1]) in ['line', 'segment', 'ray'] else f"line {out} = line({ins[0]}, {ins[1]});", "line"
    ),
    'ray': lambda ins, out, types: (f"line {out} = ray({ins[0]}, {ins[1]});", "ray"),
    'perpendicularbisector': lambda ins, out, types: (f"line {out} = bisector({ins[0]}, {ins[1]});", "line"),
    'linebisector': lambda ins, out, types: (f"line {out} = bisector({ins[0]}, {ins[1]});", "line"),
    'anglebisector': lambda ins, out, types: (
        f"line {out} = line({ins[1]}, {ins[1]} + unit({ins[0]} - {ins[1]}) + unit({ins[2]} - {ins[1]}));" if len(ins) == 3 else f"line {out} = bisector({ins[0]}, {ins[1]});", "line"
    ),
    'angularbisector': lambda ins, out, types: (
        f"line {out} = line({ins[1]}, {ins[1]} + unit({ins[0]} - {ins[1]}) + unit({ins[2]} - {ins[1]}));" if len(ins) == 3 else f"line {out} = bisector({ins[0]}, {ins[1]});", "line"
    ),
    'tangent': lambda ins, out, types: (f"line {out} = tangents({ins[1]}, {ins[0]})[0];", "line"),
    'polar': lambda ins, out, types: (f"line {out} = polar_ggb({ins[0]}, {ins[1]});", "line"),
    'orthogonalline': lambda ins, out, types: (f"line {out} = perpendicular({ins[0]}, {ins[1]});", "line"),
    'perpendicularline': lambda ins, out, types: (f"line {out} = perpendicular({ins[0]}, {ins[1]});", "line"),
}

# 3. CIRCLE COMMANDS
CIRCLE_COMMANDS = {
    'circle': lambda ins, out, types: (
        f"circle {out} = circle({ins[0]}, {ins[1]}, {ins[2]});" if len(ins) == 3 else
        (f"circle {out} = circle_point({ins[0]}, {ins[1]});" if types.get(ins[1]) == 'point' else f"circle {out} = circle({ins[0]}, {ins[1]});"),
        "circle"
    ),
    'incircle': lambda ins, out, types: (f"circle {out} = incircle(triangle({ins[0]}, {ins[1]}, {ins[2]}));", "circle"),
    'circumcircle': lambda ins, out, types: (f"circle {out} = circumcircle(triangle({ins[0]}, {ins[1]}, {ins[2]}));", "circle"),
    'circulararc': lambda ins, out, types: (f"arc {out} = arc(circle({ins[0]}, abs({ins[0]} - {ins[1]})), {ins[1]}, {ins[2]});", "arc"),
    'circumcirculararc': lambda ins, out, types: (f"arc {out} = arc_ggb({ins[0]}, {ins[1]}, {ins[2]});", "arc"),
    'circumcirclearc': lambda ins, out, types: (f"arc {out} = arc_ggb({ins[0]}, {ins[1]}, {ins[2]});", "arc"),
    'semicircle': lambda ins, out, types: (f"arc {out} = arc(circle(midpoint(segment({ins[0]}, {ins[1]})), abs({ins[0]} - {ins[1]})/2.0), {ins[1]}, {ins[0]});", "arc"),
}

# 4. CURVE, POLYGON & LOCUS COMMANDS
CURVE_COMMANDS = {
    'conic': lambda ins, out, types: (f"conic {out} = conic({ins[0]}, {ins[1]}, {ins[2]}, {ins[3]}, {ins[4]});", "conic"),
    'polygon': lambda ins, out, types: (f"path {out} = {'--'.join(ins)}--cycle;", "polygon"),
    'ellipse': lambda ins, out, types: (f"ellipse {out} = ellipse({ins[0]}, {ins[1]}, {ins[2]});", "conic"),
    'hyperbola': lambda ins, out, types: (f"hyperbola {out} = hyperbola({ins[0]}, {ins[1]}, {ins[2]});", "conic"),
    'parabola': lambda ins, out, types: (f"parabola {out} = parabola({ins[0]}, {ins[1]});", "conic"),
    
    # Standard GGB Locus equivalents mapped to 1:1 Conics
    'locuspl': lambda ins, out, types: (f"conic {out} = locus_pl({ins[0]}, {ins[1]});", "conic"),
    'locuspc': lambda ins, out, types: (f"conic {out} = locus_pc({ins[0]}, {ins[1]});", "conic"),
    'locuslc': lambda ins, out, types: (f"conic {out} = locus_lc({ins[0]}, {ins[1]});", "conic"),
    'locuscc': lambda ins, out, types: (f"conic {out} = locus_cc({ins[0]}, {ins[1]});", "conic"),
}

# 5. TRANSFORM COMMANDS
TRANSFORM_COMMANDS = {
    'mirror': lambda ins, out, types: map_transform(ins, out, types, 'mirror'),
    'translate': lambda ins, out, types: map_transform(ins, out, types, 'translate'),
    'rotate': lambda ins, out, types: map_transform(ins, out, types, 'rotate'),
    'dilate': lambda ins, out, types: map_transform(ins, out, types, 'dilate'),
}

# 6. MY CUSTOM TOOLS (User-developed tools mapped with 1:1 dynamic behavior)
MY_CUSTOM_TOOLS = {
    'harmonicconjugate3p': lambda ins, out, types: (f"point {out} = harmonic_conjugate({ins[0]}, {ins[1]}, {ins[2]});", "point"),
    'harmonicconjugateratio': lambda ins, out, types: (f"point {out} = harmonic_conjugate({ins[0]}, {ins[1]}, {ins[2]});", "point"),
    'circlecentertoline': lambda ins, out, types: (f"circle {out} = circle({ins[0]}, distance({ins[0]}, {ins[1]}));", "circle"),
    'circlecentertocircle': lambda ins, out, types: (f"circle {out} = circle_center_to_circle({ins[0]}, {ins[1]});", "circle"),
    'excircle': lambda ins, out, types: (f"circle {out} = excircle({ins[0]}, {ins[1]}, {ins[2]});", "circle"),
    'mixtilinearcircle': lambda ins, out, types: (f"circle {out} = mixtilinear_incircle({ins[0]}, {ins[1]}, {ins[2]});", "circle"),
    'exmixtilinearcircle': lambda ins, out, types: (f"circle {out} = exmixtilinear_incircle({ins[0]}, {ins[1]}, {ins[2]});", "circle"),
    'radicalline': lambda ins, out, types: (f"line {out} = radicalline({ins[0]}, {ins[1]});", "line"),
    'limitingpoints': lambda ins, out, types: (f"point {out} = limiting_point({ins[0]}, {ins[1]});", "point"),
    'homotheticcenter': lambda ins, out, types: (f"point {out} = homothetic_center({ins[0]}, {ins[1]});", "point"),
}

# Construct main mapper dictionary
COMMAND_MAPPERS = {}
for group in [POINT_COMMANDS, LINE_COMMANDS, CIRCLE_COMMANDS, CURVE_COMMANDS, TRANSFORM_COMMANDS, MY_CUSTOM_TOOLS]:
    COMMAND_MAPPERS.update(group)


class GgbToAsyConverter:
    def __init__(self, xml_content):
        self.root = ET.fromstring(xml_content)
        self.types = {}
        self.name_map = {}
        self.type_counters = defaultdict(int)

        self.part1_free_points = []
        self.part2_dependent_geom = []
        self.parsed_commands = []
        self.draw_groups = defaultdict(list)
        self.element_styles = {}
        self.label_map = {}

    def _get_new_name(self, original_label, typ):
        if typ == 'point':
            safe_name = original_label.replace("'", "_prime").replace(" ", "_")
            self.label_map[safe_name] = original_label
            return safe_name
        
        prefix_map = {
            'line': 'li',
            'segment': 'sg',
            'ray': 'ray',
            'circle': 'circ',
            'arc': 'arc',
            'conic': 'con',
            'polygon': 'poly'
        }
        prefix = prefix_map.get(typ, 'obj')
        idx = self.type_counters[prefix]
        self.type_counters[prefix] += 1
        return f"{prefix}_{idx}"

    def _get_output_labels(self, cmd_node):
        output_node = cmd_node.find('output')
        if output_node is None:
            return []
        labels = []
        i = 0
        while True:
            a = output_node.attrib.get(f'a{i}')
            if a is None:
                break
            labels.append(a)
            i += 1
        return labels

    def _generate_locus_code(self, dep_pt, driver_pt, out_name):
        dependent_vars = {driver_pt}
        subgraph_cmds = []
        
        for pcmd in getattr(self, 'parsed_commands', []):
            if any(inp in dependent_vars for inp in pcmd['ins']):
                dependent_vars.update(pcmd['out'])
                subgraph_cmds.append(pcmd['code'])

        driver_parent = "unknown_path"
        for pcmd in getattr(self, 'parsed_commands', []):
            if driver_pt in pcmd['out'] and len(pcmd['ins']) > 0:
                driver_parent = pcmd['ins'][0]
                break
                
        driver_parent_type = self.types.get(driver_parent, 'path')

        lines = []
        lines.append(f"guide generate_locus_{out_name}({driver_parent_type} driver_obj, int n=100) {{")
        lines.append(f"    point[] tmp;")
        lines.append(f"    for (int i=0; i<=n; ++i) {{")
        lines.append(f"        point {driver_pt} = relpoint(driver_obj, (real)i/n);")
        
        for c in subgraph_cmds:
            for cl in c.split('\n'):
                lines.append(f"        {cl}")
        
        lines.append(f"        tmp.push({dep_pt});")
        lines.append(f"    }}")
        lines.append(f"    guide g;")
        lines.append(f"    if (tmp.length > 0) g = tmp[0];")
        lines.append(f"    for (int i=1; i<tmp.length; ++i) g = g -- tmp[i];")
        lines.append(f"    return g;")
        lines.append(f"}}")
        lines.append(f"guide {out_name} = generate_locus_{out_name}({driver_parent});")
        
        return "\n".join(lines)

    def _map_command(self, cmd_name, original_inputs, original_outputs, point_coords):
        if not original_outputs or not any(out.strip() for out in original_outputs):
            return
            
        orig_out = original_outputs[0]
        if not orig_out.strip():
            for out in original_outputs:
                if out.strip():
                    orig_out = out
                    break
        
        # Check dependencies using the mapped input names
        for inp in original_inputs:
            mapped_inp = self.name_map.get(inp, inp)
            if mapped_inp and mapped_inp[0].isalpha() and mapped_inp not in self.types:
                cmd = f"// Skipped GGB command {cmd_name}({', '.join(original_inputs)}) -> {orig_out} due to missing dependency: {inp}"
                self.part2_dependent_geom.append(cmd)
                return


        lower_cmd = cmd_name.lower()
        
        if lower_cmd == 'locus':
            dep_pt = self.name_map.get(original_inputs[0], original_inputs[0])
            driver_pt = self.name_map.get(original_inputs[1], original_inputs[1])
            new_out = self._get_new_name(orig_out, 'guide')
            self.name_map[orig_out] = new_out
            self.types[new_out] = 'guide'
            
            locus_code = self._generate_locus_code(dep_pt, driver_pt, new_out)
            self.part2_dependent_geom.append(locus_code)
            return

        # Handle simple coordinate points
        if lower_cmd == 'point' and orig_out in point_coords:
            new_out = self._get_new_name(orig_out, 'point')
            self.name_map[orig_out] = new_out
            self.types[new_out] = 'point'
            parent = None
            if len(original_inputs) > 0:
                parent = self.name_map.get(original_inputs[0], original_inputs[0])
                parent_type = self.types.get(parent, 'path')
                suffix = 'circle' if parent_type == 'circle' else ('line' if parent_type in ['line', 'segment', 'ray', 'axis'] else 'path')
                if hasattr(self, 'coord_ratio') and self.coord_ratio != 1.0:
                    cmd = f"point {new_out} = relpoint_ggb_{suffix}({parent}, ({float(point_coords[orig_out][0]) * self.coord_ratio:.5f}, {float(point_coords[orig_out][1]) * self.coord_ratio:.5f}));"
                else:
                    cmd = f"point {new_out} = relpoint_ggb_{suffix}({parent}, ({point_coords[orig_out][0]}, {point_coords[orig_out][1]}));"
            else:
                if hasattr(self, 'coord_ratio') and self.coord_ratio != 1.0:
                    cmd = f"point {new_out} = ({float(point_coords[orig_out][0]) * self.coord_ratio:.5f}, {float(point_coords[orig_out][1]) * self.coord_ratio:.5f});"
                else:
                    cmd = f"point {new_out} = ({point_coords[orig_out][0]}, {point_coords[orig_out][1]});"
            self.part2_dependent_geom.append(cmd)
            self.parsed_commands.append({'out': [new_out], 'ins': [parent] if parent is not None else [], 'code': cmd})
            return

        # Handle dual outputs of AngularBisector
        if lower_cmd in ['angularbisector', 'anglebisector'] and len(original_inputs) == 2:
            mapped_inputs = [self.name_map.get(inp, inp) for inp in original_inputs]
            
            # Internal bisector
            out1_orig = original_outputs[0]
            out1_new = self._get_new_name(out1_orig, 'line')
            self.name_map[out1_orig] = out1_new
            self.types[out1_new] = 'line'
            self.part2_dependent_geom.append(f"line {out1_new} = bisector({mapped_inputs[0]}, {mapped_inputs[1]});")
            
            # External bisector
            if len(original_outputs) > 1:
                out2_orig = original_outputs[1]
                out2_new = self._get_new_name(out2_orig, 'line')
                self.name_map[out2_orig] = out2_new
                self.types[out2_new] = 'line'
                self.part2_dependent_geom.append(f"line {out2_new} = perpendicular(intersectionpoint({mapped_inputs[0]}, {mapped_inputs[1]}), {out1_new});")
                self.parsed_commands.append({'out': [out2_new], 'ins': mapped_inputs, 'code': f"line {out2_new} = perpendicular(intersectionpoint({mapped_inputs[0]}, {mapped_inputs[1]}), {out1_new});"})
            
            self.parsed_commands.append({'out': [out1_new], 'ins': mapped_inputs, 'code': f"line {out1_new} = bisector({mapped_inputs[0]}, {mapped_inputs[1]});"})
            return

        # Handle multiple outputs of Tangent
        if lower_cmd == 'tangent':
            mapped_inputs = [self.name_map.get(inp, inp) for inp in original_inputs]
            type_0 = self.types.get(mapped_inputs[0], 'point')
            type_1 = self.types.get(mapped_inputs[1], 'point')
            
            conic_types = ['circle', 'conic', 'ellipse', 'parabola', 'hyperbola']
            
            if type_0 == 'circle' and type_1 == 'circle':
                tgt_func = f"tangents_circle_circle({mapped_inputs[0]}, {mapped_inputs[1]})"
            elif type_0 in conic_types and type_1 in conic_types:
                tgt_func = f"tangents_ggb({mapped_inputs[0]}, {mapped_inputs[1]})"
            elif type_0 in conic_types:
                tgt_func = f"tangents({mapped_inputs[0]}, {mapped_inputs[1]})"
            else:
                tgt_func = f"tangents({mapped_inputs[1]}, {mapped_inputs[0]})"
                
            tmp_name = f"tmp_tgts_{self.type_counters.get('obj', 0)}"
            self.type_counters['obj'] = self.type_counters.get('obj', 0) + 1
            self.part2_dependent_geom.append(f"line[] {tmp_name} = {tgt_func};")
            
            for i, out_orig in enumerate(original_outputs):
                out_new = self._get_new_name(out_orig, 'line')
                self.name_map[out_orig] = out_new
                self.types[out_new] = 'line'
                self.part2_dependent_geom.append(f"line {out_new};")
                self.part2_dependent_geom.append(f"if ({tmp_name}.length > {i}) {out_new} = {tmp_name}[{i}];")
                self.parsed_commands.append({'out': [out_new], 'ins': mapped_inputs, 'code': f"line {out_new}; if ({tmp_name}.length > {i}) {out_new} = {tmp_name}[{i}];"})
            self.parsed_commands.append({'out': [tmp_name], 'ins': mapped_inputs, 'code': f"line[] {tmp_name} = {tgt_func};"})
            return

        mapped_inputs = [self.name_map.get(inp, inp) for inp in original_inputs]
        
        # Dispatch command via grouped mappers
        if lower_cmd in COMMAND_MAPPERS:
            dummy_out = "{OUT}"
            cmd_template, typ = COMMAND_MAPPERS[lower_cmd](mapped_inputs, dummy_out, self.types)
            
            new_out = self._get_new_name(orig_out, typ)
            self.name_map[orig_out] = new_out
            self.types[new_out] = typ
            
            cmd = cmd_template.replace(dummy_out, new_out)
            self.part2_dependent_geom.append(cmd)
            self.parsed_commands.append({'out': [new_out], 'ins': mapped_inputs, 'code': cmd})
        else:
            cmd = f"// Unsupported GGB command: {cmd_name}({', '.join(original_inputs)}) -> {orig_out}"
            self.part2_dependent_geom.append(cmd)

    def parse(self):
        construction = self.root.find('.//construction')
        if construction is None:
            return ""

        point_coords = {}
        for el in construction.findall('element'):
            label = el.attrib.get('label')
            
            show_node = el.find('show')
            is_visible = show_node is not None and show_node.attrib.get('object', 'false') == 'true'
            self.element_styles[label] = {'visible': is_visible}
            
            if el.attrib.get('type') == 'point':
                coords = el.find('coords')
                if coords is not None:
                    x = float(coords.attrib.get('x', '0'))
                    y = float(coords.attrib.get('y', '0'))
                    z = float(coords.attrib.get('z', '1'))
                    if abs(z) > 1e-12:
                        x /= z
                        y /= z
                    point_coords[label] = (x, y)
        
        command_outputs = set()
        for cmd in construction.findall('command'):
            command_outputs.update(self._get_output_labels(cmd))
        for expr in construction.findall('expression'):
            lbl = expr.attrib.get('label')
            if lbl:
                command_outputs.add(lbl)

        # Part 1: Free Points
        for child in construction:
            if child.tag == 'element':
                label = child.attrib.get('label')
                if label not in command_outputs:
                    new_label = self._get_new_name(label, 'point')
                    self.name_map[label] = new_label
                    self.types[new_label] = 'point'

        coord_ratio = 1.0
        if point_coords:
            xs = [float(v[0]) for v in point_coords.values()]
            ys = [float(v[1]) for v in point_coords.values()]
            w = max(xs) - min(xs) if xs else 0
            h = max(ys) - min(ys) if ys else 0
            max_dim = max(w, h)
            if max_dim > 20:
                coord_ratio = 10.0 / max_dim
        self.coord_ratio = coord_ratio

        for orig_out in self.name_map.keys():
            new_out = self.name_map[orig_out]
            typ = self.types[new_out]
            if typ == 'point' and orig_out in point_coords:
                x, y = point_coords[orig_out]
                if coord_ratio != 1.0:
                    self.part1_free_points.append(f"point {new_out} = ({float(x) * coord_ratio:.5f}, {float(y) * coord_ratio:.5f});")
                else:
                    self.part1_free_points.append(f"point {new_out} = ({x}, {y});")

        for child in construction:
            if child.tag == 'command':
                cmd_name = child.attrib.get('name')
                input_node = child.find('input')
                orig_inputs = []
                if input_node is not None:
                    i = 0
                    while True:
                        a = input_node.attrib.get(f'a{i}')
                        if a is None: break
                        orig_inputs.append(a)
                        i += 1
                        
                outputs = self._get_output_labels(child)
                self._map_command(cmd_name, orig_inputs, outputs, point_coords)
            elif child.tag == 'expression':
                import re
                label = child.attrib.get('label')
                exp = child.attrib.get('exp')
                typ = child.attrib.get('type', 'point')
                if label and exp:
                    new_label = self._get_new_name(label, typ)
                    self.name_map[label] = new_label
                    self.types[new_label] = typ
                    
                    deps = [m for m in re.findall(r'[a-zA-Z_]\w*', exp) if m in self.name_map]
                    
                    asy_exp = exp
                    for d in sorted(set(deps), key=len, reverse=True):
                        asy_exp = re.sub(rf'\b{d}\b', self.name_map[d], asy_exp)
                        
                    cmd = f"{typ} {new_label} = {asy_exp};"
                    self.part2_dependent_geom.append(cmd)
                    
                    mapped_inputs = [self.name_map[d] for d in deps]
                    self.parsed_commands.append({'out': [new_label], 'ins': mapped_inputs, 'code': cmd})

        # Part 3: Drawing & Styling
        for orig_label, style in self.element_styles.items():
            if style['visible'] and orig_label in self.name_map:
                mapped_name = self.name_map[orig_label]
                if mapped_name in self.types:
                    t = self.types[mapped_name]
                    if t == 'point':
                        self.draw_groups['point'].append(mapped_name)
                    else:
                        self.draw_groups[t].append(mapped_name)

        draw_section = ["\n// --- 3. DRAWING & STYLING ---"]
        
        # Group paths FIRST (so they are drawn underneath)
        for typ, items in self.draw_groups.items():
            if typ == 'point': continue
            if not items:
                continue
            
            if typ == 'vector':
                path_expr = " ^^ ".join([f"{item}_path" for item in items])
                pen_var = f"pen_{typ}"
                if typ not in ['line', 'segment', 'ray', 'circle', 'arc', 'conic', 'polygon', 'vector']:
                    pen_var = "black"
                draw_section.append(f"draw({path_expr}, {pen_var}, Arrow);")
            else:
                path_expr = " ^^ ".join(items)
                pen_var = f"pen_{typ}"
                if typ not in ['line', 'segment', 'ray', 'circle', 'arc', 'conic', 'polygon', 'guide']:
                    pen_var = "black"
                draw_section.append(f"draw({path_expr}, {pen_var});")

        # Point drawing phase (always drawn last so they appear on top)
        if self.draw_groups['point']:
            pts_array = " ^^ ".join(self.draw_groups['point'])
            draw_section.append(f"dot({pts_array}, Fill(white));")
            draw_section.append(f"point[] all_pts = new point[]{{{', '.join(self.draw_groups['point'])}}};")
            for pt in self.draw_groups['point']:
                orig_name = self.label_map.get(pt, pt)
                draw_section.append(f'label("${orig_name}$", {pt}, auto_align({pt}, all_pts), pen_label);')

        unitsize = 1.0
        asy_header_dynamic = ASY_HEADER.replace('unitsize(1cm);', f'unitsize({unitsize}cm);')

        # Dynamically generate styles based on what is actually drawn
        dynamic_styles = []
        if self.draw_groups['point']:
            dynamic_styles.append("pen pen_point = black;")
            dynamic_styles.append("pen pen_label = black;")
        for typ in self.draw_groups:
            if typ == 'point': continue
            color = "magenta" if typ == 'conic' else ("deepcyan" if typ in ['circle', 'arc'] else ("red" if typ == 'guide' else "black"))
            dynamic_styles.append(f"pen pen_{typ} = {color};")
            
        asy_header_dynamic = asy_header_dynamic.replace('// DYNAMIC_STYLES_HERE', '\n'.join(dynamic_styles))

        # Collect all generated code to scan for helpers
        all_generated_code = "\n".join(self.part1_free_points + self.part2_dependent_geom + draw_section)
        
        # Iteratively inject helpers to satisfy nested dependencies
        active_helpers = set()
        while True:
            added_new = False
            for helper_name, helper_code in ASY_HELPERS.items():
                if helper_name not in active_helpers:
                    used_in_main = helper_name in all_generated_code
                    used_in_helpers = any(helper_name in ASY_HELPERS[h] for h in active_helpers)
                    if used_in_main or used_in_helpers:
                        active_helpers.add(helper_name)
                        added_new = True
            if not added_new:
                break
                
        # Build the dynamic helper section in order
        dynamic_helpers_code = "\n\n".join(ASY_HELPERS[h] for h in ASY_HELPERS if h in active_helpers)
        asy_header_dynamic = asy_header_dynamic.replace('// DYNAMIC_HELPERS_HERE', dynamic_helpers_code)

        result = [asy_header_dynamic]
        result.append("// --- 1. FREE POINTS ---")
        result.extend(self.part1_free_points)
        result.append("\n// --- 2. DEPENDENT GEOMETRY ---")
        result.extend(self.part2_dependent_geom)
        result.extend(draw_section)
        
        result.append(f"\nshipout(bbox({unitsize}mm, invisible));")
        return "\n".join(result)


def convert_ggb_to_asy(ggb_path, asy_path):
    with zipfile.ZipFile(ggb_path, 'r') as z:
        xml_content = z.read('geogebra.xml').decode('utf-8')
        
    converter = GgbToAsyConverter(xml_content)
    asy_content = converter.parse()
    
    with open(asy_path, 'w', encoding='utf-8') as f:
        f.write(asy_content)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python ggb_to_asy.py <input.ggb> <output.asy>")
        sys.exit(1)
        
    in_path = sys.argv[1]
    out_path = sys.argv[2]
    try:
        convert_ggb_to_asy(in_path, out_path)
        print(f"Successfully converted {in_path} to {out_path}")
    except Exception as e:
        print(f"Error converting: {e}")
