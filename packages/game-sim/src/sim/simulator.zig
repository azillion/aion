const std = @import("std");
const types = @import("core");

pub const Simulator = struct {
    allocator: std.mem.Allocator,
    state: types.SystemState,
    time_scale: f64 = 1.0,
    next_body_id: u64 = 0,

    pub fn create(allocator: std.mem.Allocator, initial_state: types.SystemState) !*Simulator {
        const sim = try allocator.create(Simulator);
        sim.* = .{ .allocator = allocator, .state = initial_state, .time_scale = 1.0, .next_body_id = 0 };
        return sim;
    }

    pub fn destroy(self: *Simulator) void {
        // free string fields for bodies then free slice
        for (self.state.bodies) |b| {
            self.allocator.free(b.id);
            self.allocator.free(b.name);
        }
        self.allocator.free(self.state.bodies);

        // free string fields for ships then free slice
        for (self.state.ships) |s| {
            self.allocator.free(s.body.id);
            self.allocator.free(s.body.name);
        }
        self.allocator.free(self.state.ships);

        self.allocator.destroy(self);
    }
};
