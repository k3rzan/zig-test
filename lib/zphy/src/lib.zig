const std = @import("std");

pub const PolygonShapeError = error{
    PolygonIsNotSimple,
    InvalidVerticesLength,
    NotEnoughPoints,
    OutOfMemory,
    CouldNotSlicePolygon,
};

pub const Vec2 = struct {
    x: f32,
    y: f32,

    pub fn zero() Vec2 {
        return .{
            .y = 0.0,
            .x = 0.0,
        };
    }
};

pub const AABB = struct {
    min_x: f32,
    min_y: f32,
    max_x: f32,
    max_y: f32,
    is_colliding: bool,
};

pub const Collision = struct {
    bodies: []AnyRigidBody,
};

pub fn RectangleShape() type {
    return struct {
        const Self = @This();
        vertices: [4]Vec2,
        width: f32,
        height: f32,
        center: Vec2,
        is_colliding: bool,

        pub fn getCenter(v: [4]Vec2) Vec2 {
            var vector_sum: Vec2 = undefined;

            for (0..v.len) |i| {
                if (i + 1 < v.len) {
                    const current = v[i];

                    if (i == 0) {
                        vector_sum = current;
                    }

                    const next = v[i + 1];

                    vector_sum = .{
                        .x = vector_sum.x + next.x,
                        .y = vector_sum.y + next.y,
                    };
                }
            }

            return .{
                .x = vector_sum.x / v.len,
                .y = vector_sum.y / v.len,
            };
        }

        pub fn init(half_w: f32, half_h: f32) Self {
            const w = half_w * 2.0;
            const h = half_h * 2.0;

            const v = [4]Vec2{
                .{
                    .x = 0.0,
                    .y = 0.0,
                },
                .{
                    .x = 0.0,
                    .y = h,
                },
                .{
                    .x = w,
                    .y = h,
                },
                .{
                    .x = w,
                    .y = 0.0,
                },
            };
            const center = Self.getCenter(v);

            return Self{
                .width = w,
                .height = h,
                .vertices = v,
                .center = center,
                .is_colliding = false,
            };
        }
    };
}

pub const CircleShape = struct {
    radius: f32,
    center: Vec2,
    is_colliding: bool,
};

pub const Polygon = struct {
    vertices: []Vec2,
};

fn crossProduct(a: Vec2, b: Vec2, c: Vec2) f32 {
    return (((b.x - a.x) * (c.y - a.y)) - ((c.x - a.x) * (b.y - a.y)));
}

pub fn isLeft(a: Vec2, b: Vec2, c: Vec2) bool {
    return crossProduct(a, b, c) > 0;
}

fn isLeftOn(a: Vec2, b: Vec2, c: Vec2) bool {
    return crossProduct(a, b, c) >= 0;
}

fn isRight(a: Vec2, b: Vec2, c: Vec2) bool {
    return crossProduct(a, b, c) < 0;
}

pub fn isRightOn(a: Vec2, b: Vec2, c: Vec2) bool {
    return crossProduct(a, b, c) <= 0;
}

pub fn polygonAt(polygon: []Vec2, i: i32) Vec2 {
    const s: isize = @intCast(polygon.len);
    const index: isize = @intCast(i);

    if (i < 0) {
        const mod1: usize = @intCast(@rem(index, s) + s);
        return polygon[mod1];
    }
    const mod2: usize = @intCast(@rem(index, s));
    return polygon[mod2];
}

fn polygonIsReflex(polygon: []Vec2, i: i32) bool {
    const a = polygonAt(polygon, i - 1);
    const b = polygonAt(polygon, i);
    const c = polygonAt(polygon, i + 1);

    return isRight(a, b, c);
}

fn polygonCopy(polygon: []Vec2, i: usize, j: usize, targetPoly: *std.ArrayList(Vec2)) ![]Vec2 {
    if (i < j) {
        // Insert all vertices from i to j
        for (i..j) |k| {
            try targetPoly.append(polygon[k]);
        }
    } else {

        // Insert vertices 0 to j
        for (0..j) |k| {
            try targetPoly.append(polygon[k]);
        }

        // Insert vertices i to end
        for (i..polygon.len) |k| {
            try targetPoly.append(polygon[k]);
        }
    }

    return targetPoly.items;
}

fn polygonGetCutEdges(allocator: std.mem.Allocator, polygon: std.ArrayList(Vec2)) !std.ArrayList([]Vec2) {
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();
    const tmp1_alloc = arena.allocator();
    const tmp2_alloc = arena.allocator();
    const tmp_poly_alloc = arena.allocator();

    var tmp1 = std.ArrayList([]Vec2).init(tmp1_alloc);
    var tmp2 = std.ArrayList([]Vec2).init(tmp2_alloc);
    var min = std.ArrayList([]Vec2).init(allocator);

    var tmp_poly = std.ArrayList(Vec2).init(tmp_poly_alloc);

    defer tmp1.deinit();
    defer tmp2.deinit();
    defer tmp_poly.deinit();

    var nDiags: i32 = std.math.maxInt(i32);

    for (0..polygon.items.len) |i| {
        if (polygonIsReflex(polygon.items, @intCast(i))) {
            for (0..polygon.items.len) |j| {
                if (polygonCanSee(polygon, @intCast(i), @intCast(j))) {
                    const copy = try polygonCopy(polygon.items, j, i, &tmp_poly);
                    var param = std.ArrayList(Vec2).init(allocator);
                    try param.insertSlice(0, copy[0..]);

                    tmp1 = try polygonGetCutEdges(allocator, param);
                    tmp2 = try polygonGetCutEdges(allocator, param);

                    param.deinit();

                    for (tmp2.items) |tmp| {
                        try tmp1.append(tmp);
                    }

                    const length: i32 = @intCast(tmp1.items.len);
                    if (length < nDiags) {
                        min = try tmp1.clone();
                        nDiags = length;
                        var vec_group = [_]Vec2{
                            polygonAt(polygon.items, @intCast(j)),
                            polygonAt(polygon.items, @intCast(i)),
                        };
                        try min.append(&vec_group);
                    }
                }
            }
        }
    }

    return min;
}

fn polygonSlice(allocator: std.mem.Allocator, polygon: std.ArrayList(Vec2), cutEdges: std.ArrayList([]Vec2)) !std.ArrayList([]Vec2) {
    var result = std.ArrayList([]Vec2).init(allocator);
    if (cutEdges.items.len == 0) {
        try result.append(polygon.items);
        return result;
    }

    if (cutEdges.items.len == 2) {
        const polys = std.ArrayList([]Vec2).init(allocator);

        for (cutEdges.items) |cutEdge| {
            // Cut all polys
            for (polys.items, 0..) |poly, j| {
                result = try polygonSlice(allocator, poly, cutEdge);
                // Found poly! Cut and quit
                polys.orderedRemove(j);
                try polys.append(result[0]);
                try polys.append(result[1]);
                break;
            }
        }

        return polys;
    } else {

        // Was given one edge
        const i = indexOf(polygon, cutEdges[0]);
        const j = indexOf(polygon, cutEdges[1]);

        if (i != -1 and j != -1) {
            const res = std.ArrayList([]Vec2).init(allocator);
            try res.append(try polygonCopy(polygon.items, @intCast(i), @intCast(j)));
            try res.append(try polygonCopy(polygon.items, j, i));
            return res;
        } else {
            return error.CouldNotSlicePolygon.CouldNotSlicePolygon;
        }
    }
}

fn indexOf(polygon: []Vec2, vec: Vec2) i32 {
    for (0..polygon.len) |idx| {
        if (vec.x == polygon[idx].x and vec.y == polygon[idx].y) {
            return @intCast(idx);
        }
    }
    return -1;
}

// pub fn polygonDecomp(allocator: std.mem.Allocator, polygon: []Vec2) !std.ArrayList([]Vec2) {
//     var poly = std.ArrayList(Vec2).init(allocator);
//     try poly.insertSlice(0, polygon[0..]);
//
//     const edges = try polygonGetCutEdges(allocator, poly);
//     if (edges.items.len > 0) {
//         return try polygonSlice(allocator, poly, edges);
//     } else {
//         poly.deinit();
//         const default = std.ArrayList([]Vec2).init(allocator);
//         try default.append(polygon);
//         return default;
//     }
// }

fn polygonAppend(polygon: *std.ArrayList(Vec2), poly: []Vec2, from: usize, to: usize) !void {
    for (from..to) |i| {
        try polygon.append(poly[i]);
    }
}

pub fn polyDecompose(
    allocator: std.mem.Allocator,
    polygon: std.ArrayList(Vec2),
    acc_result: ?std.ArrayList([]Vec2),
    acc_reflexVertices: ?std.ArrayList(Vec2),
    acc_steinerPoints: ?std.ArrayList(Vec2),
    delta: ?u32,
    level: ?u32,
) !std.ArrayList([]Vec2) {
    const maxlevel: u32 = 100;
    var local_level: u32 = if (level) |res| res else 0;
    const dt: u32 = if (delta) |res| res else 25;
    var result = if (acc_result) |res| res else std.ArrayList([]Vec2).init(allocator);
    var reflexVertices = if (acc_reflexVertices) |res| res else std.ArrayList(Vec2).init(allocator);
    var steinerPoints = if (acc_steinerPoints) |res| res else std.ArrayList(Vec2).init(allocator);

    var upperInt: Vec2 = undefined;
    var lowerInt: Vec2 = undefined;
    var p: Vec2 = undefined;

    var upperDist: f32 = 0;
    var lowerDist: f32 = 0;
    var d: f32 = 0;
    var closestDist: f32 = 0;

    var upperIndex: usize = 0;
    var lowerIndex: usize = 0;
    var closestIndex: usize = 0;

    var lowerPoly = std.ArrayList(Vec2).init(allocator);
    var upperPoly = std.ArrayList(Vec2).init(allocator);

    const poly = try polygon.clone();

    if (poly.items.len < 3) {
        //std.debug.assert(false);
        return result;
    }

    local_level += 1;
    if (local_level > maxlevel) {
        //std.debug.assert(false);
        return result;
    }

    // std.debug.assert(poly.items.len == 5);

    for (0..poly.items.len) |index| {
        const i: i32 = @intCast(index);
        const is_reflex = polygonIsReflex(poly.items, i);

        // if (index == 0) {
        //     std.debug.assert(is_reflex == false);
        // }
        //
        // if (index == 1) {
        //     std.debug.assert(is_reflex);
        // }

        if (is_reflex) {
            try reflexVertices.append(poly.items[index]);

            upperDist = std.math.floatMax(f32);
            lowerDist = std.math.floatMax(f32);

            for (0..poly.items.len) |index_two| {
                const j: i32 = @intCast(index_two);
                if (isLeft(polygonAt(poly.items, i - 1), polygonAt(poly.items, i), polygonAt(poly.items, j)) and isRightOn(polygonAt(poly.items, i - 1), polygonAt(poly.items, i), polygonAt(poly.items, j - 1))) { // if line intersects with an edge
                    p = getIntersectionPoint(
                        polygonAt(poly.items, i - 1),
                        polygonAt(poly.items, i),
                        polygonAt(poly.items, j),
                        polygonAt(poly.items, j - 1),
                        dt,
                    ); // find the point of intersection

                    //std.debug.assert(false);
                    if (isRight(polygonAt(poly.items, i + 1), polygonAt(poly.items, i), p)) { // make sure it's inside the poly
                        d = sqdist(poly.items[index], p);
                        if (d < lowerDist) { // keep only the closest intersection
                            lowerDist = d;
                            lowerInt = p;
                            lowerIndex = index_two;
                        }
                    }
                }
                const condition: bool = isLeft(polygonAt(poly.items, i + 1), polygonAt(poly.items, i), polygonAt(poly.items, j + 1)) and isRightOn(polygonAt(poly.items, i + 1), polygonAt(poly.items, i), polygonAt(poly.items, j));

                // if (i == 1 and j == 3) {
                //     std.debug.print("there must be something wrong with the isLeft function\n", .{});
                //     std.debug.assert(condition);
                // }

                if (condition) {
                    p = getIntersectionPoint(
                        polygonAt(poly.items, i + 1),
                        polygonAt(poly.items, i),
                        polygonAt(poly.items, j),
                        polygonAt(poly.items, j + 1),
                        dt,
                    );

                    //std.debug.assert(false);

                    if (isLeft(polygonAt(poly.items, i - 1), polygonAt(poly.items, i), p)) {
                        d = sqdist(poly.items[index], p);
                        if (d < upperDist) {
                            upperDist = d;
                            upperInt = p;
                            upperIndex = index_two;
                        }
                    }
                }
            }

            const mod: usize = @rem((upperIndex + 1), poly.items.len);

            if (lowerIndex == mod) {
                // std.debug.assert(false);
                p.x = (lowerInt.x + upperInt.x) / 2;
                p.y = (lowerInt.y + upperInt.y) / 2;
                try steinerPoints.append(p);

                if (i < upperIndex) {
                    try polygonAppend(&lowerPoly, poly.items, index, upperIndex + 1);
                    try lowerPoly.append(p);
                    try upperPoly.append(p);
                    if (lowerIndex != 0) {
                        try polygonAppend(&upperPoly, poly.items, lowerIndex, poly.items.len);
                    }
                    try polygonAppend(&upperPoly, poly.items, 0, index + 1);
                } else {
                    if (index != 0) try polygonAppend(
                        &lowerPoly,
                        poly.items,
                        index,
                        poly.items.len,
                    );
                    try lowerPoly.append(p);
                    try polygonAppend(&lowerPoly, poly.items, 0, upperIndex + 1);
                    try lowerPoly.append(p);
                    try upperPoly.append(p);
                    try polygonAppend(&upperPoly, poly.items, lowerIndex, index + 1);
                }
            } else {
                std.debug.print("connect to the closest point within the triangle, i: {d}\n", .{index});

                if (lowerIndex > upperIndex) {
                    upperIndex += poly.items.len;
                }
                // std.debug.assert(lowerIndex == 6 and upperIndex == 10);
                closestDist = std.math.floatMax(f32);

                if (lowerIndex > upperIndex) {
                    //std.debug.assert(false);
                    return result;
                }

                for (lowerIndex..upperIndex + 1) |index_two| {
                    const j: i32 = @intCast(index_two);
                    if (isLeftOn(
                        polygonAt(poly.items, i - 1),
                        polygonAt(poly.items, i),
                        polygonAt(poly.items, j),
                    ) and isRightOn(
                        polygonAt(poly.items, i + 1),
                        polygonAt(poly.items, i),
                        polygonAt(poly.items, j),
                    )) {
                        d = sqdist(polygonAt(poly.items, i), polygonAt(poly.items, j));
                        if (d < closestDist and polygonCanSee(poly, i, j)) {
                            closestDist = d;
                            closestIndex = @rem(index_two, poly.items.len);
                        }
                    }
                }

                // std.debug.assert(closestDist == 13577);
                // std.debug.assert(index < closestIndex);

                if (index < closestIndex) {
                    std.debug.print("it should enter here\n", .{});
                    try polygonAppend(&lowerPoly, poly.items, index, closestIndex + 1);
                    if (closestIndex != 0) {
                        try polygonAppend(
                            &upperPoly,
                            poly.items,
                            closestIndex,
                            poly.items.len,
                        );
                    }
                    try polygonAppend(
                        &upperPoly,
                        poly.items,
                        0,
                        index + 1,
                    );
                    std.debug.print("this is lowerPoly: {any}\n", .{lowerPoly});
                    std.debug.print("this is upperPoly: {any}\n", .{upperPoly});
                } else {
                    std.debug.print("it shouldn't enter here yet\n", .{});
                    // std.debug.assert(false);
                    if (index != 0) {
                        try polygonAppend(
                            &lowerPoly,
                            poly.items,
                            index,
                            poly.items.len,
                        );
                    }
                    try polygonAppend(
                        &lowerPoly,
                        poly.items,
                        0,
                        closestIndex + 1,
                    );
                    try polygonAppend(
                        &upperPoly,
                        poly.items,
                        closestIndex,
                        index + 1,
                    );
                }
            }

            // solve smallest poly first
            if (lowerPoly.items.len < upperPoly.items.len) {
                std.debug.print("these are the lens: lower {d}, upper:{d}", .{ lowerPoly.items.len, upperPoly.items.len });
                // std.debug.assert(false);
                //std.debug.assert(false);
                const lp_res = try polyDecompose(allocator, lowerPoly, result, reflexVertices, steinerPoints, dt, local_level);
                if (lp_res.items.len > 0) result = lp_res;
                const up_res = try polyDecompose(allocator, upperPoly, result, reflexVertices, steinerPoints, dt, local_level);
                if (up_res.items.len > 0) result = up_res;
            } else {
                std.debug.print("Entered correctly this is result so far: {any}!\n", .{result});
                std.debug.print("And this is upperPoly so far: {any}\n", .{upperPoly});
                // std.debug.assert(result.items.len == 0);
                // std.debug.assert(false);
                //std.debug.assert(false);
                const up_res = try polyDecompose(allocator, upperPoly, result, reflexVertices, steinerPoints, dt, local_level);
                if (up_res.items.len > 0) result = up_res;
                const lp_res = try polyDecompose(allocator, lowerPoly, result, reflexVertices, steinerPoints, dt, local_level);
                if (lp_res.items.len > 0) result = lp_res;
            }
            //std.debug.assert(false);
            return result;
        }
    }

    // std.debug.assert(poly.items.len > 0);
    try result.append(poly.items);
    //std.debug.assert(false);
    std.debug.print("this is result so far: {any}\n", .{result.items});
    // std.debug.assert(false);
    return result;
}

fn safe_index_substract(index: usize, num: usize) usize {
    return if (index == 0) 0 else index - num;
}

fn scalar_eq(a: f32, b: f32, precision: u32) bool {
    const p: f32 = @floatFromInt(precision);
    return @abs(a - b) <= p;
}

fn getIntersectionPoint(p1: Vec2, p2: Vec2, q1: Vec2, q2: Vec2, dt: u32) Vec2 {
    const delta = dt;
    const a1: f32 = p2.y - p1.y;
    const b1: f32 = p1.x - p2.x;
    const c1: f32 = (a1 * p1.x) + (b1 * p1.y);
    const a2: f32 = q2.y - q1.y;
    const b2: f32 = q1.x - q2.x;
    const c2: f32 = (a2 * q1.x) + (b2 * q1.y);
    const det: f32 = (a1 * b2) - (a2 * b1);

    if (!scalar_eq(det, 0, delta)) {
        return .{ .x = ((b2 * c1) - (b1 * c2)) / det, .y = ((a1 * c2) - (a2 * c1)) / det };
    } else {
        return .{ .x = 0.0, .y = 0.0 };
    }
}

fn lineInt(l1: [2]Vec2, l2: [2]Vec2) Vec2 {
    var i = Vec2{ .x = 0.0, .y = 0.0 };

    const a1 = l1[1].y - l1[0].y;
    const b1 = l1[0].x - l1[1].x;

    const c1 = a1 * l1[0].x + b1 * l1[0].y;

    const a2 = l2[1].y - l2[0].y;
    const b2 = l2[0].x - l2[1].x;

    const c2 = a2 * l2[0].x + b2 * l2[0].y;

    const det = a1 * b2 - a2 * b1;

    if (!(@abs(det) == 0)) { // lines are not parallel
        i.x = (b2 * c1 - b1 * c2) / det;
        i.y = (a1 * c2 - a2 * c1) / det;
    }
    return i;
}

fn sqdist(a: Vec2, b: Vec2) f32 {
    const dx = b.x - a.x;
    const dy = b.y - a.y;
    return dx * dx + dy * dy;
}

fn polygonCanSee(polygon: std.ArrayList(Vec2), a: i32, b: i32) bool {
    var l1: [2]Vec2 = undefined;
    var l2: [2]Vec2 = undefined;

    if (isLeftOn(polygonAt(polygon.items, a + 1), polygonAt(polygon.items, a), polygonAt(polygon.items, b)) and isRightOn(polygonAt(polygon.items, a - 1), polygonAt(polygon.items, a), polygonAt(polygon.items, b))) {
        return false;
    }
    const dist = sqdist(polygonAt(polygon.items, a), polygonAt(polygon.items, b));

    for (0..polygon.items.len) |i| { // for each edge
        var index_a: usize = 0;
        const index: i32 = @intCast(i);
        if (a > 0) index_a = @intCast(a);
        if (@rem((i + 1), polygon.items.len) == index_a or i == index_a) { // ignore incident edges
            continue;
        }
        if (isLeftOn(polygonAt(polygon.items, a), polygonAt(polygon.items, b), polygonAt(polygon.items, index + 1)) and isRightOn(polygonAt(polygon.items, a), polygonAt(polygon.items, b), polygonAt(polygon.items, index))) { // if diag intersects an edge
            l1 = [2]Vec2{ polygonAt(polygon.items, a), polygonAt(polygon.items, b) };
            l2 = [2]Vec2{ polygonAt(polygon.items, index), polygonAt(polygon.items, index + 1) };

            const p = lineInt(l1, l2);

            if (sqdist(polygonAt(polygon.items, a), p) < dist) { // if edge is blocking visibility to b
                return false;
            }
        }
    }

    return true;
}

fn lineSegmentsIntersect(p1: Vec2, p2: Vec2, q1: Vec2, q2: Vec2) bool {
    const dx = p2.x - p1.x;
    const dy = p2.y - p1.y;
    const da = q2.x - q1.x;
    const db = q2.y - q1.y;

    // segments are parallel
    if ((da * dy - db * dx) == 0) {
        return false;
    }

    const s = (dx * (q1.y - p1.y) + dy * (p1.x - q1.x)) / (da * dy - db * dx);
    const t = (da * (p1.y - q1.y) + db * (q1.x - p1.x)) / (db * dx - da * dy);

    return (s >= 0 and s <= 1 and t >= 0 and t <= 1);
}

pub fn polygonReverse(allocator: std.mem.Allocator, polygon: std.ArrayList(Vec2)) ![]Vec2 {
    var tmp = std.ArrayList(Vec2).init(allocator);
    var poly = try polygon.clone();
    const N = poly.items.len;

    for (0..N) |_| {
        try tmp.append(poly.pop());
    }

    for (0..N) |i| {
        try poly.insert(i, tmp.items[i]);
    }

    return poly.items;
}

pub fn polygonMakeCCW(allocator: std.mem.Allocator, polygon: std.ArrayList(Vec2)) ![]Vec2 {
    var br: i32 = 0;

    // find bottom right point
    for (polygon.items, 0..) |v, i| {
        const new_br: usize = if (br > 0) @intCast(br) else 0;
        if (v.y < polygon.items[new_br].y or (v.y == polygon.items[new_br].y and v.x > polygon.items[new_br].x)) {
            br = @intCast(i);
        }
    }

    // reverse poly if clockwise
    if (!isLeft(polygonAt(polygon.items, br - 1), polygonAt(polygon.items, br), polygonAt(polygon.items, br + 1))) {
        return try polygonReverse(allocator, polygon);
    } else {
        return polygon.items;
    }
}

pub fn isPolygonSimple(polygon: []Vec2) bool {
    // Check
    var ilen = if (polygon.len > 1) polygon.len - 2 else 0;
    for (0..ilen) |i| {
        const j_len = if (i > 0) i - 1 else 0;
        for (0..j_len) |j| {
            if (lineSegmentsIntersect(polygon[i], polygon[i + 1], polygon[j], polygon[j + 1])) {
                return false;
            }
        }
    }

    // Check the segment between the last and the first point to all others
    ilen = if (polygon.len > 1) polygon.len - 2 else 0;
    const len = if (polygon.len > 0) polygon.len - 1 else 0;
    // //std.debug.assert(false);
    for (1..ilen) |i| {
        if (lineSegmentsIntersect(polygon[0], polygon[len], polygon[i], polygon[i + 1])) {
            return false;
        }
    }

    return true;
}

const PolygonShape = struct {
    const Self = @This();
    vertices: [][]Vec2,

    pub fn createFromVertices(allocator: std.mem.Allocator, vertices: []Vec2) !Self {
        if (vertices.len < 3) return PolygonShapeError.InvalidVerticesLength;
        var polygon = std.ArrayList(Vec2).init(allocator);
        defer polygon.deinit();
        try polygon.insertSlice(0, vertices[0..]);

        if (isPolygonSimple(polygon.items)) {
            const ccw_p = try polygonMakeCCW(allocator, polygon);

            var ccw_polygon = std.ArrayList(Vec2).init(allocator);
            try ccw_polygon.insertSlice(0, ccw_p);

            const sliced_polygon = try polyDecompose(allocator, ccw_polygon, null, null, null, null, null);

            return .{
                .vertices = sliced_polygon.items,
            };
        }

        return PolygonShapeError.PolygonIsNotSimple;
    }
};

const ShapeUnion = union(enum) {
    circle: CircleShape,
    rectangle: RectangleShape(),
    polygon: PolygonShape,
};

pub fn RigidBody(comptime s: type) type {
    return struct {
        const Self = @This();

        shape: s,
        position: Vec2,
        aabb: AABB,
        velocity: Vec2,
        is_colliding: bool,

        pub fn init(
            shape: s,
            position: Vec2,
        ) Self {
            var aabb: AABB = undefined;
            switch (s) {
                CircleShape => {
                    const circle: CircleShape = shape;
                    aabb = AABB{
                        .min_x = (position.x + circle.center.x - circle.radius),
                        .min_y = (position.y + circle.center.y - circle.radius),
                        .max_x = (position.x + circle.center.x + circle.radius),
                        .max_y = (position.y + circle.center.y + circle.radius),
                        .is_colliding = false,
                    };
                },
                RectangleShape() => {
                    const rectangle: RectangleShape() = shape;

                    aabb = AABB{
                        .min_x = position.x,
                        .min_y = position.y,
                        .max_x = (position.x + rectangle.width),
                        .max_y = (position.y + rectangle.height),
                        .is_colliding = false,
                    };
                },
                else => @panic("shape type not supported"),
            }
            std.debug.print("Body initialized correctly\n", .{});

            return .{
                .shape = shape,
                .position = position,
                .aabb = aabb,
                .velocity = Vec2.zero(),
                .is_colliding = false,
            };
        }

        pub fn setVelocity(self: *Self, v: Vec2) void {
            self.velocity = v;
        }
    };
}

const WorldOptions = struct {
    enable_debug_draw: bool,
};

const AnyRigidBody = union(enum) {
    circle: RigidBody(CircleShape),
    rectangle: RigidBody(RectangleShape()),

    fn getActiveType(self: @This()) type {
        return switch (self) {
            inline else => |b| @TypeOf(b),
        };
    }

    fn getActiveShapeCenter(self: @This()) Vec2 {
        return switch (self) {
            inline else => |b| b.shape.center,
        };
    }

    fn getActiveShapeVertices(self: @This()) ?[4]Vec2 {
        return switch (self) {
            .circle => return null,
            .rectangle => |b| return b.shape.vertices,
        };
    }

    fn getActiveShape(self: @This()) type {
        return switch (self) {
            inline else => |b| @TypeOf(b.shape),
        };
    }

    fn getActiveAABB(self: @This()) AABB {
        return switch (self) {
            inline else => |b| b.aabb,
        };
    }

    fn getActiveShapeValue(self: @This()) getActiveShape() {
        return switch (self) {
            inline else => |b| b.shape,
        };
    }
};

fn getCollisionDisplacement(v1: Vec2, v2: Vec2, v3: Vec2, v4: Vec2) ?Vec2 {
    const line_a = Vec2{ .x = v2.x - v1.x, .y = v2.y - v1.y };
    const line_b = Vec2{ .x = v4.x - v3.x, .y = v4.y - v3.y };

    const denom = (line_a.x * line_b.y) - (line_a.y * line_b.x);
    const u: f32 = ((v3.x - v1.x) * line_a.y - (v3.y - v1.y) * line_a.x) / denom;
    const t: f32 = ((v3.x - v1.x) * line_b.y - (v3.y - v1.y) * line_b.x) / denom;

    if (0.0 <= u and u <= 1.0 and 0.0 <= t and t <= 1.0) {
        return .{
            .x = (1.0 - t) * (v1.x - v2.x),
            .y = (1.0 - t) * (v1.y - v2.y),
        };
    }

    return null;
}

pub fn World() type {
    return struct {
        const Self = @This();

        r_bodies: std.ArrayList(AnyRigidBody),

        draw_circle_shape: ?*const fn (radius: f32, center: Vec2) void,
        draw_rectangle_shape: ?*const fn (vertices: [4]Vec2, center: Vec2, is_colliding: bool) void,
        draw_polygon_shape: ?*const fn (polygon: [][]Vec2, position: Vec2) void,
        draw_aabb: ?*const fn (aabb: AABB) void,
        allocator: std.mem.Allocator,

        pub fn init(allocator: std.mem.Allocator) Self {
            return .{
                .allocator = allocator,
                .r_bodies = std.ArrayList(AnyRigidBody).init(allocator),
                .draw_circle_shape = null,
                .draw_rectangle_shape = null,
                .draw_polygon_shape = null,
                .draw_aabb = null,
            };
        }

        pub fn addRigidBody(self: *Self, body: AnyRigidBody) !void {
            try self.r_bodies.append(body);
        }

        pub fn step(self: *Self, delta_time: f32) void {
            self.update_bodies(delta_time);
            self.resolve_collisions(delta_time);
        }

        fn resolve_collisions(self: *Self, _: f32) void {
            var collisions = std.ArrayList(Collision).init(self.allocator);

            //this isn't efficient, but it's what I have for now
            //later search for "spatial partitioning"
            for (self.r_bodies.items, 0..) |receiver, idx| {
                var buffer = std.ArrayList(AnyRigidBody).init(self.allocator);
                const receiver_aabb = receiver.getActiveAABB();

                var is_colliding: bool = false;

                for (self.r_bodies.items, 0..) |collider, jdx| {
                    if (jdx == idx) continue;

                    const collider_aabb = collider.getActiveAABB();

                    if (receiver_aabb.max_x >= collider_aabb.min_x and receiver_aabb.min_y <= collider_aabb.max_y and receiver_aabb.max_y >= collider_aabb.min_y and receiver_aabb.min_x <= collider_aabb.max_x) {
                        std.debug.print("\n", .{});
                        std.debug.print("AABBs TOUCHING now idx:({d}) jdx: ({d})\n", .{ idx, jdx });
                        switch (collider) {
                            .circle => self.r_bodies.items[jdx].circle.aabb.is_colliding = true,
                            .rectangle => self.r_bodies.items[jdx].rectangle.aabb.is_colliding = true,
                        }
                        is_colliding = true;
                        buffer.append(collider) catch std.debug.print("there's an error\n", .{});

                        const receiver_pos = receiver.rectangle.position;
                        const collider_pos = collider.rectangle.position;

                        const receiver_center = receiver.getActiveShapeCenter();
                        const receiver_vertices = receiver.getActiveShapeVertices();
                        const collider_vertices = receiver.getActiveShapeVertices();

                        var is_intersected: bool = false;

                        if (receiver_vertices) |r_vertices| {
                            if (collider_vertices) |c_vertices| {
                                for (r_vertices, 0..) |rv, i| {
                                    if (is_intersected) break;

                                    const next_rv = r_vertices[@rem(i + 1, r_vertices.len)];

                                    for (c_vertices, 0..) |cv, j| {
                                        const index = @rem(j + 1, c_vertices.len);
                                        const next_cv = c_vertices[index];
                                        //
                                        const middle_rv = Vec2{
                                            .x = (rv.x + next_rv.x) / 2.0,
                                            .y = (rv.y + next_rv.y) / 2.0,
                                        };

                                        const v1: Vec2 = .{
                                            .x = receiver_pos.x + receiver_center.x,
                                            .y = receiver_pos.y + receiver_center.y,
                                        };
                                        const v2: Vec2 = .{
                                            .x = receiver_pos.x + middle_rv.x,
                                            .y = receiver_pos.y + middle_rv.y,
                                        };
                                        const v3: Vec2 = .{
                                            .x = collider_pos.x + cv.x,
                                            .y = collider_pos.y + cv.y,
                                        };

                                        const v4: Vec2 = .{
                                            .x = collider_pos.x + next_cv.x,
                                            .y = collider_pos.y + next_cv.y,
                                        };

                                        if (getCollisionDisplacement(v1, v2, v3, v4)) |displacement| {
                                            is_intersected = true;

                                            std.debug.print("is intersecting\n", .{});
                                            std.debug.print("this is the displacement: ({d},{d})\n", .{ displacement.x, displacement.y });
                                            switch (receiver) {
                                                .circle => self.r_bodies.items[idx].circle.shape.is_colliding = true,
                                                .rectangle => {
                                                    std.debug.print("substracting position to the receiver\n", .{});
                                                    self.r_bodies.items[idx].rectangle.shape.is_colliding = true;
                                                },
                                            }
                                            switch (collider) {
                                                .circle => self.r_bodies.items[jdx].circle.shape.is_colliding = true,
                                                .rectangle => {
                                                    self.r_bodies.items[jdx].rectangle.shape.is_colliding = true;
                                                    self.r_bodies.items[jdx].rectangle.position.x -= displacement.x;
                                                    self.r_bodies.items[jdx].rectangle.position.y -= displacement.y;
                                                    self.r_bodies.items[jdx].rectangle.aabb.min_x -= displacement.x;
                                                    self.r_bodies.items[jdx].rectangle.aabb.min_y -= displacement.y;
                                                    self.r_bodies.items[jdx].rectangle.aabb.max_x -= displacement.x;
                                                    self.r_bodies.items[jdx].rectangle.aabb.max_y -= displacement.y;
                                                },
                                            }
                                            break;
                                        } else {
                                            is_intersected = false;
                                            switch (receiver) {
                                                .circle => self.r_bodies.items[idx].circle.shape.is_colliding = false,
                                                .rectangle => self.r_bodies.items[idx].rectangle.shape.is_colliding = false,
                                            }
                                            switch (collider) {
                                                .circle => self.r_bodies.items[jdx].circle.shape.is_colliding = false,
                                                .rectangle => self.r_bodies.items[jdx].rectangle.shape.is_colliding = false,
                                            }
                                            std.debug.print("is NOT intersecting\n", .{});
                                        }
                                    }
                                }
                            }
                        }
                    } else {
                        std.debug.print("\n", .{});
                        std.debug.print("AABBs NOT touching idx:({d}) jdx: ({d})\n", .{ idx, jdx });
                        switch (collider) {
                            .circle => self.r_bodies.items[jdx].circle.aabb.is_colliding = false,
                            .rectangle => {
                                self.r_bodies.items[jdx].rectangle.aabb.is_colliding = false;
                                self.r_bodies.items[jdx].rectangle.shape.is_colliding = false;
                            },
                        }

                        switch (receiver) {
                            .circle => self.r_bodies.items[idx].circle.aabb.is_colliding = false,
                            .rectangle => {
                                self.r_bodies.items[idx].rectangle.aabb.is_colliding = false;
                                self.r_bodies.items[idx].rectangle.shape.is_colliding = false;
                            },
                        }
                    }
                }

                // if (is_colliding) {
                //     switch (receiver) {
                //         .circle => self.r_bodies.items[idx].circle.aabb.is_colliding = true,
                //         .rectangle => self.r_bodies.items[idx].rectangle.aabb.is_colliding = true,
                //     }
                //     buffer.append(receiver) catch std.debug.print("there's an error\n", .{});
                // }

                if (buffer.items.len > 0) {
                    collisions.append(.{ .bodies = buffer.items }) catch std.debug.print("another error\n", .{});
                }
            }

            // if (collisions.items.len > 0) {
            //     std.debug.print("these are the collisions: {any}\n", .{collisions.items});
            // } else {
            //     std.debug.print("There are no collisions\n", .{});
            // }

            // for (collisions.items) |collision| {
            //     for (collision.bodies) |body| {
            //         switch (body) {
            //             .rectangle => |b| {
            //                 b.shape.vertices
            //             }
            //         }
            //     }
            // }
        }

        fn update_bodies(self: *Self, delta_time: f32) void {
            for (self.r_bodies.items, 0..) |item, idx| {
                switch (item) {
                    .rectangle => |*body| {
                        const current_pos = self.r_bodies.items[idx].rectangle.position;
                        const current_aabb = self.r_bodies.items[idx].rectangle.aabb;

                        self.r_bodies.items[idx].rectangle.position.x = current_pos.x + body.velocity.x * delta_time;
                        self.r_bodies.items[idx].rectangle.position.y = current_pos.y + body.velocity.y * delta_time;
                        self.r_bodies.items[idx].rectangle.aabb = .{
                            .min_x = current_aabb.min_x + body.velocity.x * delta_time,
                            .min_y = current_aabb.min_y + body.velocity.y * delta_time,
                            .max_x = current_aabb.max_x + body.velocity.x * delta_time,
                            .max_y = current_aabb.max_y + body.velocity.y * delta_time,
                            .is_colliding = current_aabb.is_colliding,
                        };
                    },
                    .circle => |*body| {
                        const current_pos = self.r_bodies.items[idx].circle.position;
                        const current_aabb = self.r_bodies.items[idx].circle.aabb;

                        self.r_bodies.items[idx].circle.position.x = current_pos.x + body.velocity.x * delta_time;
                        self.r_bodies.items[idx].circle.position.y = current_pos.y + body.velocity.y * delta_time;

                        self.r_bodies.items[idx].circle.aabb = .{
                            .min_x = current_aabb.min_x + body.velocity.x * delta_time,
                            .min_y = current_aabb.min_y + body.velocity.y * delta_time,
                            .max_x = current_aabb.max_x + body.velocity.x * delta_time,
                            .max_y = current_aabb.max_y + body.velocity.y * delta_time,
                            .is_colliding = false,
                        };
                    },
                }
            }
        }

        pub fn draw(self: *Self) void {
            for (self.r_bodies.items) |body| {
                switch (body) {
                    .circle => |val| {
                        if (self.draw_circle_shape) |callback| {
                            callback(val.shape.radius, .{
                                .y = val.position.y + val.shape.center.y,
                                .x = val.position.x + val.shape.center.x,
                            });
                        }

                        if (self.draw_aabb) |callback| {
                            callback(val.aabb);
                        }
                    },
                    .rectangle => |val| {
                        if (self.draw_rectangle_shape) |callback| {
                            var shape_vertices: [4]Vec2 = [4]Vec2{
                                Vec2{
                                    .x = 0.0,
                                    .y = 0.0,
                                },
                                Vec2{
                                    .x = 0.0,
                                    .y = 0.0,
                                },
                                Vec2{
                                    .x = 0.0,
                                    .y = 0.0,
                                },
                                Vec2{
                                    .x = 0.0,
                                    .y = 0.0,
                                },
                            };

                            for (val.shape.vertices, 0..) |v, i| {
                                shape_vertices[i] = .{
                                    .x = val.position.x + v.x,
                                    .y = val.position.y + v.y,
                                };
                            }

                            callback(shape_vertices, .{
                                .y = val.position.y + val.shape.center.y,
                                .x = val.position.x + val.shape.center.x,
                            }, val.shape.is_colliding);
                        }

                        if (self.draw_aabb) |callback| {
                            callback(val.aabb);
                        }
                    },
                }
            }
        }
    };
}
