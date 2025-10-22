import type { Authority, InputState } from './authority';
import { quat, vec3 } from 'gl-matrix';
import type { Body, SystemState, Vec3, Ship } from '../shared/types';
import { G } from '../shared/constants';

export class LocalAuthority implements Authority {
	private state: SystemState;
	private timeScalar: number;

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
          };

		const playerShip: Ship = {
			id: 'player-ship',
			name: 'AION-1',
			// Place ship in a 10,000km altitude orbit around Earth
			// Position is Earth's position plus an offset for the orbit
			position: [
				earth.position[0] + (earth.radius + 10000),
				earth.position[1],
				earth.position[2]
			],
			// Velocity is Earth's velocity plus the tangential orbital velocity
			velocity: [
				earth.velocity[0],
				earth.velocity[1] + Math.sqrt(G * earth.mass / (earth.radius + 10000)),
				earth.velocity[2]
			],
			radius: 0.1,
			mass: 1e6,
			albedo: [0.8, 0.8, 0.9],
			emissive: [0.1, 0.3, 1.0],
			orientation: [0, 0, 0, 1],
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

		// Update player ship orientation from input deltas (stable FPS-style: world-yaw, local-pitch)
		if ((input.deltaX !== 0 || input.deltaY !== 0)) {
			const ship = this.state.bodies.find(b => b.id === 'player-ship') as Ship | undefined;
			if (ship) {
				const sensitivity = 0.0015; // radians per pixel
				const yawAngle = -input.deltaX * sensitivity;
				const pitchAngle = -input.deltaY * sensitivity;
				const q = quat.fromValues(ship.orientation[0], ship.orientation[1], ship.orientation[2], ship.orientation[3]);
				// 1) World yaw: rotate around world Y (pre-multiply)
				if (yawAngle !== 0) {
					const yawQ = quat.create();
					quat.setAxisAngle(yawQ, vec3.fromValues(0, 1, 0), yawAngle);
					quat.multiply(q, yawQ, q);
				}
				// 2) Local pitch: rotate around camera right (post-multiply)
				if (pitchAngle !== 0) {
					const right = vec3.create();
					vec3.transformQuat(right, vec3.fromValues(1, 0, 0), q);
					const pitchQ = quat.create();
					quat.setAxisAngle(pitchQ, right, pitchAngle);
					quat.multiply(q, q, pitchQ);
				}
				quat.normalize(q, q);
				ship.orientation = [q[0], q[1], q[2], q[3]];
			}
		}

		// Thrust & maneuvering
		{
			const ship = this.state.bodies.find(b => b.id === 'player-ship') as Ship | undefined;
			if (ship) {
				// Instant face Sun (KeyF): rotate -Z forward to direction toward sun
				if (input.keys.has('KeyF')) {
					const sun = this.state.bodies.find(b => b.id === 'sol');
					if (sun) {
						const toSun = vec3.fromValues(
							sun.position[0] - ship.position[0],
							sun.position[1] - ship.position[1],
							sun.position[2] - ship.position[2]
						);
						vec3.normalize(toSun, toSun);
						const baseForward = vec3.fromValues(0, 0, -1);
						const qrot = quat.create();
						quat.rotationTo(qrot, baseForward, toSun);
						ship.orientation = [qrot[0], qrot[1], qrot[2], qrot[3]];
					}
				}

				const THRUST_FORCE = 5e5; // N (0.5 m/s^2 for 1e6 kg)
				const MANEUVER_FORCE = THRUST_FORCE / 4;
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
				if (input.keys.has('KeyS')) { vec3.scaleAndAdd(totalForceVec, totalForceVec, shipForward, -THRUST_FORCE); }
				if (input.keys.has('KeyD')) { vec3.scaleAndAdd(totalForceVec, totalForceVec, shipRight, MANEUVER_FORCE); }
				if (input.keys.has('KeyA')) { vec3.scaleAndAdd(totalForceVec, totalForceVec, shipRight, -MANEUVER_FORCE); }
				if (input.keys.has('Space')) { vec3.scaleAndAdd(totalForceVec, totalForceVec, shipUp, MANEUVER_FORCE); }
				if (input.keys.has('ControlLeft')) { vec3.scaleAndAdd(totalForceVec, totalForceVec, shipUp, -MANEUVER_FORCE); }

				// Braking (KeyX)
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
	}

	public setTimeScale(scale: number): void {
		this.timeScalar = scale;
	}

	public addBody(body: Omit<Body, 'id'>): void {
		const uniqueId = `${Date.now().toString(36)}-${Math.random().toString(36).slice(2, 10)}`;
		const fullBody: Body = { id: uniqueId, ...body };
		this.state.bodies.push(fullBody);
	}
}