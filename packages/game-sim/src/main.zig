// This file is the root of the game-sim module.
// It exposes all public submodules.
pub const core = @import("core");
pub const math = @import("math");
pub const physics = @import("physics");
pub const sim = @import("sim");
pub const ffi = @import("ffi");

// Empty main, as this is a library module for WASM.
pub fn main() void {}
