const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});

    const optimize = b.standardOptimizeOption(.{});

    const module = b.addModule("zphy", .{
        .root_source_file = b.path("src/lib.zig"),
        .optimize = optimize,
        .target = target,
    });

    const main_lib = b.addExecutable(.{
        .name = "main",
        .target = target,
        .optimize = optimize,
        .root_source_file = b.path("src/main.zig"),
    });

    b.installArtifact(main_lib);

    main_lib.root_module.addImport("zphy", module);

    const run_cmd = b.addRunArtifact(main_lib);

    run_cmd.step.dependOn(b.getInstallStep());
    // const zbox_step = zbox.builder.getInstallStep();
    // run_cmd.step.dependOn(zbox_step);

    if (b.args) |args| {
        run_cmd.addArgs(args);
    }

    const run_step = b.step("run", "Run the app");
    run_step.dependOn(&run_cmd.step);

    const run = b.addTest(.{
        .root_source_file = b.path("src/main.zig"),
        .target = target,
    });

    const test_step = b.step("test", "Run the test file");
    const tests = b.addTest(.{
        .root_source_file = b.path("src/test.zig"),
        .target = target,
    });

    tests.root_module.addImport("zphy", module);
    run.root_module.addImport("zphy", module);

    const run_unit_tests = b.addRunArtifact(tests);
    test_step.dependOn(&run_unit_tests.step);
}
