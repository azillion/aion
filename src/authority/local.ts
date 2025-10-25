import type { Authority, InputState } from './authority';
import { quat, vec3 } from 'gl-matrix';
import type { Body, SystemState, Vec3, Ship, TerrainParams } from '../shared/types';
import { G } from '../shared/constants';

export class LocalAuthority implements Authority {
	private state: SystemState;
	private timeScalar: number;
	private autoLanding: { active: boolean; targetId: string | null } = { active: false, targetId: null };

	constructor() {
		const sun: Body = {
			id: 'sol',
			name: 'Sun',
			position: [0, 0, 0],
			velocity: [0, 0, 0],
			radius: 696340,
			mass: 1.989e30,
			albedo: [1.0, 1.0, 0.9],
			emissive: [1000.0, 1000.0, 900.0],
		};

		const earth_dist = 149.6e6; // 1 AU in km
		const earth_vel = Math.sqrt(G * sun.mass / earth_dist); // Circular orbit velocity

		const earth: Body = {
			id: 'earth',
			name: 'Earth',
			// Start Earth on the Y axis to give the camera a side-on, lit view
			position: [0, earth_dist, 0],
			// Velocity tangential to position for circular orbit
			velocity: [-earth_vel, 0, 0],
			radius: 6371,
			mass: 5.972e24,
			albedo: [0.2, 0.3, 0.8],
			terrain: { radius: 6371, seaLevel: 1.5, maxHeight: 0.0013, noiseSeed: 42.0 } as TerrainParams,
		};
		
		const moon_dist = 384400; // Lunar distance from Earth
		const moon_vel = Math.sqrt(G * earth.mass / moon_dist);

        const moon: Body = {
            id: 'moon',
            name: 'Moon',
            position: [0, earth_dist - moon_dist, 0],
            velocity: [-earth_vel + moon_vel, 0, 0],
            radius: 1737,
            mass: 7.347e22,
            albedo: [0.5, 0.5, 0.5],
            terrain: { radius: 1737, seaLevel: 0.0, maxHeight: 0.0015, noiseSeed: 1337.0 } as TerrainParams,
          };

		const playerShip: Ship = {
			id: 'player-ship',
			name: 'AION-1',
			// Place ship in a ~35,786 km (GEO) altitude circular orbit around Earth
			// Position is Earth's position plus an offset for the orbit
			position: [
				earth.position[0] + (earth.radius + 35786),
				earth.position[1],
				earth.position[2]
			],
			// Velocity is Earth's velocity plus the tangential orbital velocity
			velocity: [
				earth.velocity[0],
				earth.velocity[1] + Math.sqrt(G * earth.mass / (earth.radius + 35786)),
				earth.velocity[2]
			],
			radius: 0.1,
			mass: 1e6,
			albedo: [0.8, 0.8, 0.9],
			emissive: [0.1, 0.3, 1.0],
			orientation: [0, 0, 0, 1],
			angularVelocity: [0, 0, 0],
			thrust: [0, 0, 0],
		};

		this.state = {
			timestamp: Date.now(),
			bodies: [sun, earth, moon, playerShip],
		};

		// Default to 1.0 for real-time simulation (1 simulated second per 1 real second)
		this.timeScalar = 1.0;
	}

	async query(): Promise<SystemState> {
		return this.state;
	}

	async tick(deltaTime: number, input: InputState): Promise<void> {
		const dt = deltaTime * this.timeScalar;

		// Direct mouse-look deprecated: orientation now controlled by RCS torques

		// Reset frame flags
		this.state.flags = { precision: false, killRotation: false };

		// Thrust & maneuvering
		{
			const ship = this.state.bodies.find(b => b.id === 'player-ship') as Ship | undefined;
			if (ship) {
				// Prograde Lock (KeyP): face velocity vector relative to Sun (temporary frame)
				if (input.keys.has('KeyP')) {
					const shipVel = vec3.fromValues(ship.velocity[0], ship.velocity[1], ship.velocity[2]);
					const sun = this.state.bodies.find(b => b.id === 'sol');
					if (sun) {
						const sunVel = vec3.fromValues(sun.velocity[0], sun.velocity[1], sun.velocity[2]);
						const relativeVel = vec3.subtract(vec3.create(), shipVel, sunVel);
						if (vec3.length(relativeVel) > 1) {
							vec3.normalize(relativeVel, relativeVel);
							const baseForward = vec3.fromValues(0, 0, -1);
							const qrot = quat.create();
							quat.rotationTo(qrot, baseForward, relativeVel);
							ship.orientation = [qrot[0], qrot[1], qrot[2], qrot[3]];
						}
					}
				}

				// Barycenter Lock (KeyB): face toward system barycenter (excluding the ship)
				if (input.keys.has('KeyB')) {
					let totalMass = 0;
					let cx = 0, cy = 0, cz = 0;
					for (const body of this.state.bodies) {
						if (body.id === 'player-ship') continue;
						totalMass += body.mass;
						cx += body.mass * body.position[0];
						cy += body.mass * body.position[1];
						cz += body.mass * body.position[2];
					}
					if (totalMass > 0) {
						const bx = cx / totalMass;
						const by = cy / totalMass;
						const bz = cz / totalMass;
						const toBary = vec3.fromValues(bx - ship.position[0], by - ship.position[1], bz - ship.position[2]);
						if (vec3.length(toBary) > 1) {
							vec3.normalize(toBary, toBary);
							const baseForward = vec3.fromValues(0, 0, -1);
							const qrot = quat.create();
							quat.rotationTo(qrot, baseForward, toBary);
							ship.orientation = [qrot[0], qrot[1], qrot[2], qrot[3]];
						}
					}
				}

				// Precision Mode (Tab): scale forces/torques for fine control
				const isPrecisionMode = input.keys.has('Backquote');
				if (isPrecisionMode) this.state.flags!.precision = true;
				const precisionMultiplier = isPrecisionMode ? 0.1 : 1.0;

				const THRUST_FORCE = 5e5 * precisionMultiplier; // Newtons
				const MANEUVER_FORCE = THRUST_FORCE / 4; // Maneuver thrusters scaled too
				const BASE_FORWARD = vec3.fromValues(0, 0, -1);
				const BASE_RIGHT = vec3.fromValues(1, 0, 0);
				const BASE_UP = vec3.fromValues(0, 1, 0);
				const q = quat.fromValues(ship.orientation[0], ship.orientation[1], ship.orientation[2], ship.orientation[3]);
				const shipForward = vec3.create();
				const shipRight = vec3.create();
				const shipUp = vec3.create();
				vec3.transformQuat(shipForward, BASE_FORWARD, q);
				vec3.transformQuat(shipRight, BASE_RIGHT, q);
				vec3.transformQuat(shipUp, BASE_UP, q);

				const totalForceVec = vec3.create();
				if (input.keys.has('KeyW')) { vec3.scaleAndAdd(totalForceVec, totalForceVec, shipForward, THRUST_FORCE); }
				if (input.keys.has('KeyS')) { vec3.scaleAndAdd(totalForceVec, totalForceVec, shipForward, -THRUST_FORCE / 2); }
				if (input.keys.has('KeyD')) { vec3.scaleAndAdd(totalForceVec, totalForceVec, shipRight, MANEUVER_FORCE); }
				if (input.keys.has('KeyA')) { vec3.scaleAndAdd(totalForceVec, totalForceVec, shipRight, -MANEUVER_FORCE); }
				if (input.keys.has('Space')) { vec3.scaleAndAdd(totalForceVec, totalForceVec, shipUp, MANEUVER_FORCE); }
				if (input.keys.has('ControlLeft')) { vec3.scaleAndAdd(totalForceVec, totalForceVec, shipUp, -MANEUVER_FORCE); }

				// RCS rotation: map keys to angular acceleration; Backspace kills rotation
				const ANGULAR_ACCELERATION = 0.5 * precisionMultiplier; // rad/s^2
				const angVel = vec3.fromValues(ship.angularVelocity[0], ship.angularVelocity[1], ship.angularVelocity[2]);
				if (input.keys.has('Backspace')) { // Rotational dampening not affected by precision mode
					this.state.flags!.killRotation = true;
					const DAMPENING_FACTOR = 0.95;
					vec3.scale(angVel, angVel, DAMPENING_FACTOR);
					if (vec3.length(angVel) < 0.001) {
						angVel[0] = 0; angVel[1] = 0; angVel[2] = 0;
					}
				} else {
					const angularTorque = vec3.create();
					// Pitch (R/F)
					if (input.keys.has('KeyR')) { angularTorque[0] -= ANGULAR_ACCELERATION; }
					if (input.keys.has('KeyF')) { angularTorque[0] += ANGULAR_ACCELERATION; }
					// Yaw (Z/C)
					if (input.keys.has('KeyZ')) { angularTorque[1] += ANGULAR_ACCELERATION; }
					if (input.keys.has('KeyC')) { angularTorque[1] -= ANGULAR_ACCELERATION; }
					// Roll (Q/E)
					if (input.keys.has('KeyQ')) { angularTorque[2] += ANGULAR_ACCELERATION; }
					if (input.keys.has('KeyE')) { angularTorque[2] -= ANGULAR_ACCELERATION; }
					// Integrate angular velocity
					vec3.scaleAndAdd(angVel, angVel, angularTorque, dt);
				}
				ship.angularVelocity = [angVel[0], angVel[1], angVel[2]];

				// Apply angular velocity to orientation (convert to degrees for fromEuler)
				const deltaQuat = quat.create();
				const degX = angVel[0] * dt * (180 / Math.PI);
				const degY = angVel[1] * dt * (180 / Math.PI);
				const degZ = angVel[2] * dt * (180 / Math.PI);
				quat.fromEuler(deltaQuat, degX, degY, degZ);
				quat.multiply(q, q, deltaQuat);
				quat.normalize(q, q);
				ship.orientation = [q[0], q[1], q[2], q[3]];

				// Braking (KeyX): counter forward/backward component
				if (input.keys.has('KeyX')) {
					const vel = vec3.fromValues(ship.velocity[0], ship.velocity[1], ship.velocity[2]);
					const vForward = vec3.dot(vel, shipForward);
					if (Math.abs(vForward) > 0.1) {
						const brakeDir = vec3.create();
						vec3.scale(brakeDir, shipForward, -Math.sign(vForward));
						vec3.scaleAndAdd(totalForceVec, totalForceVec, brakeDir, THRUST_FORCE);
					}
				}

				// Apply acceleration using scaled dt to match gravity step
				// Store thrust acceleration for substep integration
				const acc = vec3.create();
				vec3.scale(acc, totalForceVec, 1 / ship.mass);
				(ship as any)._thrustAcc = acc;
			}
		}

		// Substep integration for stability
		let remaining = dt;
		const maxStep = 10; // seconds per substep
		const playerIndex = this.state.bodies.findIndex(b => b.id === 'player-ship');
		while (remaining > 0) {
			const h = Math.min(remaining, maxStep);
			const accelerations: Vec3[] = this.state.bodies.map(() => [0, 0, 0]);
			// Gravity accelerations
			for (let i = 0; i < this.state.bodies.length; i++) {
				for (let j = i + 1; j < this.state.bodies.length; j++) {
					const bodyA = this.state.bodies[i];
					const bodyB = this.state.bodies[j];
					const dx = bodyB.position[0] - bodyA.position[0];
					const dy = bodyB.position[1] - bodyA.position[1];
					const dz = bodyB.position[2] - bodyA.position[2];
					const distSq = dx * dx + dy * dy + dz * dz;
					const dist = Math.sqrt(Math.max(distSq, 1e-9));
					const forceMagnitude = (G * bodyA.mass * bodyB.mass) / Math.max(distSq, 1e-9);
					const forceX = forceMagnitude * (dx / dist);
					const forceY = forceMagnitude * (dy / dist);
					const forceZ = forceMagnitude * (dz / dist);
					accelerations[i][0] += forceX / bodyA.mass;
					accelerations[i][1] += forceY / bodyA.mass;
					accelerations[i][2] += forceZ / bodyA.mass;
					accelerations[j][0] -= forceX / bodyB.mass;
					accelerations[j][1] -= forceY / bodyB.mass;
					accelerations[j][2] -= forceZ / bodyB.mass;
				}
			}
			// Add thrust acceleration if present
			if (playerIndex >= 0) {
				const acc: any = (this.state.bodies[playerIndex] as any)._thrustAcc;
				if (acc) {
					accelerations[playerIndex][0] += acc[0];
					accelerations[playerIndex][1] += acc[1];
					accelerations[playerIndex][2] += acc[2];
				}

				// Auto-landing guidance: apply corrective acceleration towards surface normal and damp velocity
				const ship = this.state.bodies[playerIndex] as Ship;
				if (this.autoLanding.active && this.autoLanding.targetId) {
					const target = this.state.bodies.find(b => b.id === this.autoLanding.targetId);
					if (target) {
						// Vector from planet center to ship
						const rx = ship.position[0] - target.position[0];
						const ry = ship.position[1] - target.position[1];
						const rz = ship.position[2] - target.position[2];
						const r = Math.sqrt(rx*rx + ry*ry + rz*rz);
						const nx = rx / Math.max(r, 1e-6);
						const ny = ry / Math.max(r, 1e-6);
						const nz = rz / Math.max(r, 1e-6);

						// Estimate safe surface radius: use terrain if present, else radius
						const baseR = target.terrain ? target.terrain.radius : target.radius;
						const maxH = target.terrain ? target.terrain.maxHeight * baseR : 0.0;
						const safety = Math.max(0.05 * baseR, 2.0); // 5% radius or 2 km buffer
						const desiredR = baseR + maxH + safety;

						// Radial and tangential velocity components
						const vx = ship.velocity[0] - target.velocity[0];
						const vy = ship.velocity[1] - target.velocity[1];
						const vz = ship.velocity[2] - target.velocity[2];
						const vRad = vx*nx + vy*ny + vz*nz;
						const vtx = vx - vRad*nx;
						const vty = vy - vRad*ny;
						const vtz = vz - vRad*nz;

						// Guidance: target zero tangential, small downward vRad while above surface, zero near surface
						const kTangential = 1e-4; // weak lateral damping
						accelerations[playerIndex][0] += -kTangential * vtx;
						accelerations[playerIndex][1] += -kTangential * vty;
						accelerations[playerIndex][2] += -kTangential * vtz;

						// Radial control: PD towards desiredR
						const dr = r - desiredR;
						const kPos = 5e-6;
						const kVel = 5e-4;
						const aRad = -kPos * dr - kVel * vRad; // accelerate inward when above, brake descent near surface
						accelerations[playerIndex][0] += aRad * nx;
						accelerations[playerIndex][1] += aRad * ny;
						accelerations[playerIndex][2] += aRad * nz;

						// Completion check: close and slow
						const speed = Math.sqrt(vx*vx + vy*vy + vz*vz);
						if (r <= desiredR + 0.5 && speed < 0.01) {
							this.autoLanding.active = false;
						}
					}
				}
			}
			// Integrate
			for (let i = 0; i < this.state.bodies.length; i++) {
				const body = this.state.bodies[i];
				body.velocity[0] += accelerations[i][0] * h;
				body.velocity[1] += accelerations[i][1] * h;
				body.velocity[2] += accelerations[i][2] * h;
				body.position[0] += body.velocity[0] * h;
				body.position[1] += body.velocity[1] * h;
				body.position[2] += body.velocity[2] * h;
			}
			this.state.timestamp += h * 1000;
			remaining -= h;
		}

		// Penetration prevention: clamp ship above surface after integration
		if (playerIndex >= 0) {
			const ship = this.state.bodies[playerIndex] as Ship;
			// Find nearest massive body (planet/moon) by distance
			let nearest: Body | null = null;
			let nearestR = Infinity;
			for (const b of this.state.bodies) {
				if (b.id === ship.id) continue;
				const dx = ship.position[0] - b.position[0];
				const dy = ship.position[1] - b.position[1];
				const dz = ship.position[2] - b.position[2];
				const r = Math.sqrt(dx*dx + dy*dy + dz*dz);
				if (r < nearestR) { nearestR = r; nearest = b; }
			}
			if (nearest) {
				const baseR = nearest.terrain ? nearest.terrain.radius : nearest.radius;
				const maxH = nearest.terrain ? nearest.terrain.maxHeight * baseR : 0.0;
				const safety = Math.max(0.05 * baseR, 2.0);
				const minR = baseR + maxH + safety;
				if (nearestR < minR) {
					// Push ship out to surface along radial
					const rx = ship.position[0] - nearest.position[0];
					const ry = ship.position[1] - nearest.position[1];
					const rz = ship.position[2] - nearest.position[2];
					const inv = 1.0 / Math.max(nearestR, 1e-6);
					const nx = rx * inv, ny = ry * inv, nz = rz * inv;
					ship.position[0] = nearest.position[0] + nx * minR;
					ship.position[1] = nearest.position[1] + ny * minR;
					ship.position[2] = nearest.position[2] + nz * minR;
					// Zero radial velocity into the ground
					const vx = ship.velocity[0] - nearest.velocity[0];
					const vy = ship.velocity[1] - nearest.velocity[1];
					const vz = ship.velocity[2] - nearest.velocity[2];
					const vRad = vx*nx + vy*ny + vz*nz;
					if (vRad < 0) {
						ship.velocity[0] -= vRad * nx;
						ship.velocity[1] -= vRad * ny;
						ship.velocity[2] -= vRad * nz;
					}
				}
			}
		}
	}

	public autoLand(targetBodyId: string | null): void {
		if (!targetBodyId) return;
		const target = this.state.bodies.find(b => b.id === targetBodyId);
		const ship = this.state.bodies.find(b => b.id === 'player-ship');
		if (!target || !ship) return;
		this.autoLanding.active = true;
		this.autoLanding.targetId = targetBodyId;
	}

	public setTimeScale(scale: number): void {
		this.timeScalar = Math.max(0, scale); // Ensure time doesn't go backwards
	}

	public addBody(body: Omit<Body, 'id'>): void {
		const uniqueId = `${Date.now().toString(36)}-${Math.random().toString(36).slice(2, 10)}`;
		const fullBody: Body = { id: uniqueId, ...body };
		this.state.bodies.push(fullBody);
	}
}