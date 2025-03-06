const std = @import("std");
const zphy = @import("zphy");

fn areListsEqual(expected_result: [][]zphy.Vec2, result: [][]zphy.Vec2) bool {
    if (expected_result.len != result.len) return false;

    const ilen = if (expected_result.len > 0) expected_result.len - 1 else 0;

    for (0..ilen) |idx| {
        const ier = expected_result[idx];
        const ir = result[idx];
        const jlen = if (ier.len > 0) ier.len - 1 else 0;

        for (0..jlen) |jdx| {
            std.debug.print("iteration [{d}][{d}]\n", .{ idx, jdx });
            std.debug.print("expected: [{d},{d}] vs actual: [{d},{d}]", .{
                ier[jdx].x,
                ier[jdx].y,
                ir[idx].x,
                ir[idx].y,
            });
            const are_items_equal = ier[jdx].x == ir[jdx].x and ier[jdx].y == ir[jdx].y;
            if (!are_items_equal) {
                return false;
            }
        }
    }
    return true;
}

fn getCase3Vertices(allocator: std.mem.Allocator) !std.ArrayList(zphy.Vec2) {
    var vertices = std.ArrayList(zphy.Vec2).init(allocator);

    try vertices.append(.{ .x = 160, .y = 56 });
    try vertices.append(.{ .x = 177, .y = 84 });
    try vertices.append(.{ .x = 196, .y = 61 });
    try vertices.append(.{ .x = 214, .y = 105 });
    try vertices.append(.{ .x = 316, .y = 56 });
    try vertices.append(.{ .x = 318, .y = 257 });
    try vertices.append(.{ .x = 219, .y = 255 });
    try vertices.append(.{ .x = 202, .y = 231 });
    try vertices.append(.{ .x = 184, .y = 253 });
    try vertices.append(.{ .x = 166, .y = 200 });
    try vertices.append(.{ .x = 53, .y = 252 });
    try vertices.append(.{ .x = 56, .y = 78 });

    return vertices;
}

test "simple test" {
    const rb = zphy.RigidBody(zphy.CircleShape).init(.{
        .radius = 1.0,
        .center = zphy.Vec2{ .x = 0.0, .y = 0.0 },
    }, .{ .x = 1.0, .y = 2.0 });

    try std.testing.expect(rb.shape.radius == 1.0);
    try std.testing.expect(rb.position.x == 1.0);
}

test "Check if polygon is simple" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const alloc = arena.allocator();
    var vertices = std.ArrayList(zphy.Vec2).init(alloc);

    try vertices.append(.{ .x = 0.0909091, .y = 0.818182 });
    try vertices.append(.{ .x = 5.09091, .y = 0.818182 });
    try vertices.append(.{ .x = 7.09091, .y = 1.81818 });
    try vertices.append(.{ .x = 9.09091, .y = 3.81818 });
    try vertices.append(.{ .x = 10.0909, .y = 5.81818 });
    try vertices.append(.{ .x = 10.0909, .y = 10.8182 });
    try vertices.append(.{ .x = 9.09091, .y = 12.8182 });
    try vertices.append(.{ .x = 7.09091, .y = 14.8182 });
    try vertices.append(.{ .x = 4.09091, .y = 16.8182 });
    try vertices.append(.{ .x = 4.09091, .y = 17.8182 });
    try vertices.append(.{ .x = 8.09091, .y = 19.8182 });
    try vertices.append(.{ .x = 9.09091, .y = 20.8182 });
    try vertices.append(.{ .x = 10.0909, .y = 21.8182 });
    try vertices.append(.{ .x = 10.0909, .y = 22.8182 });
    try vertices.append(.{ .x = 7.09091, .y = 22.8182 });
    try vertices.append(.{ .x = 7.09091, .y = 21.8182 });
    try vertices.append(.{ .x = -1.90909, .y = 21.8182 });
    try vertices.append(.{ .x = -1.90909, .y = 22.8182 });
    try vertices.append(.{ .x = -4.90909, .y = 22.8182 });
    try vertices.append(.{ .x = -4.90909, .y = 21.8182 });
    try vertices.append(.{ .x = -3.90909, .y = 20.8182 });
    try vertices.append(.{ .x = -2.90909, .y = 19.8182 });
    try vertices.append(.{ .x = 1.09091, .y = 17.8182 });
    try vertices.append(.{ .x = 1.09091, .y = 16.8182 });
    try vertices.append(.{ .x = -1.90909, .y = 14.8182 });
    try vertices.append(.{ .x = -3.90909, .y = 12.8182 });
    try vertices.append(.{ .x = -4.90909, .y = 10.8182 });
    try vertices.append(.{ .x = -4.90909, .y = 5.81818 });
    try vertices.append(.{ .x = -3.90909, .y = 3.81818 });
    try vertices.append(.{ .x = -1.90909, .y = 1.81818 });

    const condition = zphy.isPolygonSimple(vertices.items);
    try std.testing.expect(condition);
}

// test "reversing polygon" {
//     var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
//     defer arena.deinit();
//     const alloc = arena.allocator();
//     var vertices = std.ArrayList(zphy.Vec2).init(alloc);
//
//     try vertices.append(.{ .x = 0.0909091, .y = 0.818182 });
//     try vertices.append(.{ .x = 5.09091, .y = 0.818182 });
//     try vertices.append(.{ .x = 7.09091, .y = 1.81818 });
//     try vertices.append(.{ .x = 9.09091, .y = 3.81818 });
//     try vertices.append(.{ .x = 10.0909, .y = 5.81818 });
//     try vertices.append(.{ .x = 10.0909, .y = 10.8182 });
//     try vertices.append(.{ .x = 9.09091, .y = 12.8182 });
//     try vertices.append(.{ .x = 7.09091, .y = 14.8182 });
//     try vertices.append(.{ .x = 4.09091, .y = 16.8182 });
//     try vertices.append(.{ .x = 4.09091, .y = 17.8182 });
//     try vertices.append(.{ .x = 8.09091, .y = 19.8182 });
//     try vertices.append(.{ .x = 9.09091, .y = 20.8182 });
//     try vertices.append(.{ .x = 10.0909, .y = 21.8182 });
//     try vertices.append(.{ .x = 10.0909, .y = 22.8182 });
//     try vertices.append(.{ .x = 7.09091, .y = 22.8182 });
//     try vertices.append(.{ .x = 7.09091, .y = 21.8182 });
//     try vertices.append(.{ .x = -1.90909, .y = 21.8182 });
//     try vertices.append(.{ .x = -1.90909, .y = 22.8182 });
//     try vertices.append(.{ .x = -4.90909, .y = 22.8182 });
//     try vertices.append(.{ .x = -4.90909, .y = 21.8182 });
//     try vertices.append(.{ .x = -3.90909, .y = 20.8182 });
//     try vertices.append(.{ .x = -2.90909, .y = 19.8182 });
//     try vertices.append(.{ .x = 1.09091, .y = 17.8182 });
//     try vertices.append(.{ .x = 1.09091, .y = 16.8182 });
//     try vertices.append(.{ .x = -1.90909, .y = 14.8182 });
//     try vertices.append(.{ .x = -3.90909, .y = 12.8182 });
//     try vertices.append(.{ .x = -4.90909, .y = 10.8182 });
//     try vertices.append(.{ .x = -4.90909, .y = 5.81818 });
//     try vertices.append(.{ .x = -3.90909, .y = 3.81818 });
//     try vertices.append(.{ .x = -1.90909, .y = 1.81818 });
//
//     const polygon_shape = try zphy.PolygonShape.createFromVertices(alloc, vertices.items);
//
//     std.debug.print("these are the new: {any}\n", .{polygon_shape});
//
//     try std.testing.expect(true);
// }

test "testing isLeft and PolygonAt" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const alloc = arena.allocator();

    const vertices = try getCase3Vertices(alloc);

    const poly = try zphy.polygonMakeCCW(alloc, vertices);

    const is_left = zphy.isLeft(zphy.polygonAt(poly, 2), zphy.polygonAt(poly, 1), zphy.polygonAt(poly, 11));
    try std.testing.expect(is_left);
}

test "testing isRightOn" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const alloc = arena.allocator();

    const vertices = try getCase3Vertices(alloc);

    const poly = try zphy.polygonMakeCCW(alloc, vertices);
    const is_right_on = zphy.isRightOn(zphy.polygonAt(poly, 2), zphy.polygonAt(poly, 1), zphy.polygonAt(poly, 10));

    try std.testing.expect(is_right_on);
}

test "complex decomposition" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const alloc = arena.allocator();

    const vertices = try getCase3Vertices(alloc);

    const result = try zphy.PolygonShape.createFromVertices(alloc, vertices.items);

    std.debug.print("this is the result: {any}\n", .{result.vertices});

    var expected_result = std.ArrayList([]zphy.Vec2).init(alloc);

    try expected_result.append(@constCast(&[_]zphy.Vec2{
        .{ .x = 166.0, .y = 200.0 },
        .{ .x = 53.0, .y = 252.0 },
        .{ .x = 56.0, .y = 78.0 },
        .{ .x = 160.0, .y = 56.0 },
        .{ .x = 177.0, .y = 84.0 },
    }));
    try expected_result.append(@constCast(&[_]zphy.Vec2{
        .{ .x = 166.0, .y = 200.0 },
        .{ .x = 177.0, .y = 84.0 },
        .{ .x = 196.0, .y = 61.0 },
        .{ .x = 214.0, .y = 105.0 },
    }));
    try expected_result.append(@constCast(&[_]zphy.Vec2{
        .{ .x = 202.0, .y = 231.0 },
        .{ .x = 184.0, .y = 253.0 },
        .{ .x = 166.0, .y = 200.0 },
        .{ .x = 214.0, .y = 105.0 },
    }));
    try expected_result.append(@constCast(&[_]zphy.Vec2{
        .{ .x = 214.0, .y = 105.0 },
        .{ .x = 316.0, .y = 56.0 },
        .{ .x = 318.0, .y = 257.0 },
        .{ .x = 219.0, .y = 255.0 },
        .{ .x = 202.0, .y = 231.0 },
    }));

    const are_equal = areListsEqual(expected_result.items, result.vertices);

    try std.testing.expect(are_equal);
}
// test "simple decomposition" {
//     var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
//     defer arena.deinit();
//     const alloc = arena.allocator();
//     var vertices = std.ArrayList(zphy.Vec2).init(alloc);
//
//     try vertices.append(.{ .x = 192.0, .y = 40.0 });
//     try vertices.append(.{ .x = 32.0, .y = 48.0 });
//     try vertices.append(.{ .x = 78.0, .y = 154.0 });
//     try vertices.append(.{ .x = 118.0, .y = 108.0 });
//     try vertices.append(.{ .x = 160.0, .y = 150.0 });
//
//     const result = try zphy.PolygonShape.createFromVertices(alloc, vertices.items);
//
//     const expected_result = [2][4]zphy.Vec2{
//         .{
//             .{ .x = 116.60822510822511, .y = 43.769588744588745 },
//             .{ .x = 192, .y = 40 },
//             .{ .x = 160, .y = 150 },
//             .{ .x = 118, .y = 108 },
//         },
//         .{
//             .{ .x = 118, .y = 108 },
//             .{ .x = 78, .y = 154 },
//             .{ .x = 32, .y = 48 },
//             .{ .x = 116.60822510822511, .y = 43.769588744588745 },
//         },
//     };
//
//     try std.testing.expect(result.vertices.len == expected_result.len);
//
//     const ilen = if (expected_result.len > 0) expected_result.len - 1 else 0;
//     for (0..ilen) |idx| {
//         const ier = expected_result[idx];
//         const ir = result.vertices[idx];
//         const jlen = if (ier.len > 0) ier.len - 1 else 0;
//
//         for (0..jlen) |jdx| {
//             std.debug.print("iteration [{d}][{d}]\n", .{ idx, jdx });
//             std.debug.print("expected: [{d},{d}] vs actual: [{d},{d}]", .{
//                 ier[jdx].x,
//                 ier[jdx].y,
//                 ir[idx].x,
//                 ir[idx].y,
//             });
//             try std.testing.expect(ier[jdx].x == ir[jdx].x and ier[jdx].y == ir[jdx].y);
//         }
//     }
// }

// test "test polygon decomposition" {
//     var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
//     defer arena.deinit();
//     const alloc = arena.allocator();
//     var vertices = std.ArrayList(zphy.Vec2).init(alloc);
//
//     try vertices.append(.{ .x = 0.0909091, .y = 0.818182 });
//     try vertices.append(.{ .x = 5.09091, .y = 0.818182 });
//     try vertices.append(.{ .x = 7.09091, .y = 1.81818 });
//     try vertices.append(.{ .x = 9.09091, .y = 3.81818 });
//     try vertices.append(.{ .x = 10.0909, .y = 5.81818 });
//     try vertices.append(.{ .x = 10.0909, .y = 10.8182 });
//     try vertices.append(.{ .x = 9.09091, .y = 12.8182 });
//     try vertices.append(.{ .x = 7.09091, .y = 14.8182 });
//     try vertices.append(.{ .x = 4.09091, .y = 16.8182 });
//     try vertices.append(.{ .x = 4.09091, .y = 17.8182 });
//     try vertices.append(.{ .x = 8.09091, .y = 19.8182 });
//     try vertices.append(.{ .x = 9.09091, .y = 20.8182 });
//     try vertices.append(.{ .x = 10.0909, .y = 21.8182 });
//     try vertices.append(.{ .x = 10.0909, .y = 22.8182 });
//     try vertices.append(.{ .x = 7.09091, .y = 22.8182 });
//     try vertices.append(.{ .x = 7.09091, .y = 21.8182 });
//     try vertices.append(.{ .x = -1.90909, .y = 21.8182 });
//     try vertices.append(.{ .x = -1.90909, .y = 22.8182 });
//     try vertices.append(.{ .x = -4.90909, .y = 22.8182 });
//     try vertices.append(.{ .x = -4.90909, .y = 21.8182 });
//     try vertices.append(.{ .x = -3.90909, .y = 20.8182 });
//     try vertices.append(.{ .x = -2.90909, .y = 19.8182 });
//     try vertices.append(.{ .x = 1.09091, .y = 17.8182 });
//     try vertices.append(.{ .x = 1.09091, .y = 16.8182 });
//     try vertices.append(.{ .x = -1.90909, .y = 14.8182 });
//     try vertices.append(.{ .x = -3.90909, .y = 12.8182 });
//     try vertices.append(.{ .x = -4.90909, .y = 10.8182 });
//     try vertices.append(.{ .x = -4.90909, .y = 5.81818 });
//     try vertices.append(.{ .x = -3.90909, .y = 3.81818 });
//     try vertices.append(.{ .x = -1.90909, .y = 1.81818 });
//
//     const result = try zphy.polyDecompose(alloc, vertices, null, null, null, null, null);
//
//     std.debug.print("this is the result of the decomposition: {any}\n", .{result.items});
//
//     try std.testing.expect(false);
// }
