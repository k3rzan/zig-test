const std = @import("std");
const rl = @import("raylib");
const zphy = @import("zphy");

fn draw_circle(radius: f32, center: zphy.Vec2) void {
    rl.drawCircleLinesV(.{ .x = center.x * 50.0, .y = center.y * 50.0 }, radius * 50.0, rl.Color.red);
}

fn draw_rectangle(vertices: [4]zphy.Vec2, center: zphy.Vec2, _: bool) void {
    rl.drawCircleLinesV(.{ .x = center.x * 50.0, .y = center.y * 50.0 }, 5.0, rl.Color.sky_blue);

    for (0..vertices.len) |i| {
        const v = vertices[i];
        const vi = vertices[@rem(i + 1, vertices.len)];

        rl.drawLineV(.{ .x = v.x * 50.0, .y = v.y * 50.0 }, .{ .x = vi.x * 50.0, .y = vi.y * 50.0 }, rl.Color.sky_blue);
    }
}
fn draw_aabb(aabb: zphy.AABB) void {
    var color: rl.Color = rl.Color.sky_blue;
    if (aabb.is_colliding) {
        color = rl.Color.red;
    }
    rl.drawRectangleLinesEx(
        .{
            .x = aabb.min_x * 50.0,
            .y = aabb.min_y * 50.0,
            .width = (aabb.max_x - aabb.min_x) * 50.0,
            .height = (aabb.max_y - aabb.min_y) * 50.0,
        },
        2.0,
        color,
    );
}
pub fn main() !void {
    const screenWidth = 920;
    const screenHeight = 780;

    rl.initWindow(screenWidth, screenHeight, "raylib [core] example - basic window");
    defer rl.closeWindow();

    rl.setTargetFPS(120);

    const rectangle = zphy.RectangleShape().init(1.0, 2.0);

    var sample_body = zphy.RigidBody(zphy.RectangleShape()).init(rectangle, .{ .x = 1.0, .y = 1.0 });
    var top_body = zphy.RigidBody(zphy.RectangleShape()).init(rectangle, .{ .x = 4.0, .y = -5.0 });
    const rectangle_body = zphy.RigidBody(zphy.RectangleShape()).init(rectangle, .{ .x = 4.0, .y = 1.0 });

    var another_body = zphy.RigidBody(zphy.RectangleShape()).init(rectangle, .{ .x = 8.0, .y = 5.0 });
    const another_collider_body = zphy.RigidBody(zphy.RectangleShape()).init(rectangle, .{ .x = 2.0, .y = 5.0 });

    another_body.setVelocity(.{ .x = -1.0, .y = 0.0 });
    sample_body.setVelocity(.{ .x = 1.0, .y = 0.0 });
    top_body.setVelocity(.{ .x = 0.0, .y = 1.0 });

    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    const alloc = arena.allocator();
    defer arena.deinit();

    var world = zphy.World().init(alloc);
    // world.draw_circle_shape = draw_circle;
    // world.draw_aabb = draw_aabb;
    world.draw_rectangle_shape = draw_rectangle;

    try world.addRigidBody(.{ .rectangle = another_body });
    try world.addRigidBody(.{ .rectangle = another_collider_body });
    try world.addRigidBody(.{ .rectangle = sample_body });
    try world.addRigidBody(.{ .rectangle = top_body });
    try world.addRigidBody(.{ .rectangle = rectangle_body });

    const camera = rl.Camera2D{
        .zoom = 0.5,
        .offset = .{ .x = 0.0, .y = 0.0 },
        .target = .{ .y = -500.0, .x = -500.0 },
        .rotation = 0.0,
    };

    while (!rl.windowShouldClose()) // Detect window close button or ESC key
    {
        rl.beginDrawing();
        defer rl.endDrawing();

        rl.beginMode2D(camera);

        rl.clearBackground(rl.Color.white);

        world.draw();
        world.step(1.0 / 60.0);
    }
}
