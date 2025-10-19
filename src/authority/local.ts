import type { Authority } from './authority';
import type { Body, SystemState, Vec3 } from '../shared/types';
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
			emissive: [1.0, 1.0, 0.9],
		};

		const earth_dist = 149.6e6; // 1 AU in km
		const earth_vel = Math.sqrt(G * sun.mass / earth_dist); // Circular orbit velocity

		const earth: Body = {
			id: 'earth',
			name: 'Earth',
			position: [earth_dist, 0, 0],
			velocity: [0, earth_vel, 0], // Correct orbital velocity
			radius: 6371,
			mass: 5.972e24,
			albedo: [0.2, 0.3, 0.8],
		};
		
		const moon_dist = 384400; // Lunar distance from Earth
		const moon_vel = Math.sqrt(G * earth.mass / moon_dist);

		const moon: Body = {
			id: 'moon',
			name: 'Moon',
			position: [earth_dist + moon_dist, 0, 0],
			velocity: [0, earth_vel + moon_vel, 0], // Earth's velocity + Moon's orbital velocity
			radius: 1737,
			mass: 7.347e22,
			albedo: [0.5, 0.5, 0.5],
		};

		this.state = {
			timestamp: Date.now(),
			bodies: [sun, earth, moon],
		};

		// Default to 30 days of simulated time per real second
		this.timeScalar = 30 * 24 * 60 * 60;
	}

	async query(): Promise<SystemState> {
		return this.state;
	}

	async tick(deltaTime: number): Promise<void> {
		const dt = deltaTime * this.timeScalar;

		const accelerations: Vec3[] = this.state.bodies.map(() => [0, 0, 0]);

		// Phase 1: Calculate forces and accelerations
		for (let i = 0; i < this.state.bodies.length; i++) {
			for (let j = i + 1; j < this.state.bodies.length; j++) {
				const bodyA = this.state.bodies[i];
				const bodyB = this.state.bodies[j];

				const dx = bodyB.position[0] - bodyA.position[0];
				const dy = bodyB.position[1] - bodyA.position[1];
				const dz = bodyB.position[2] - bodyA.position[2];

				const distSq = dx * dx + dy * dy + dz * dz;
				const dist = Math.sqrt(distSq);

				const forceMagnitude = (G * bodyA.mass * bodyB.mass) / distSq;

				const forceX = forceMagnitude * (dx / dist);
				const forceY = forceMagnitude * (dy / dist);
				const forceZ = forceMagnitude * (dz / dist);

				// Apply force to body A
				accelerations[i][0] += forceX / bodyA.mass;
				accelerations[i][1] += forceY / bodyA.mass;
				accelerations[i][2] += forceZ / bodyA.mass;

				// Apply opposite force to body B
				accelerations[j][0] -= forceX / bodyB.mass;
				accelerations[j][1] -= forceY / bodyB.mass;
				accelerations[j][2] -= forceZ / bodyB.mass;
			}
		}

		// Phase 2: Update velocities and positions
		for (let i = 0; i < this.state.bodies.length; i++) {
			const body = this.state.bodies[i];
			
			body.velocity[0] += accelerations[i][0] * dt;
			body.velocity[1] += accelerations[i][1] * dt;
			body.velocity[2] += accelerations[i][2] * dt;

			body.position[0] += body.velocity[0] * dt;
			body.position[1] += body.velocity[1] * dt;
			body.position[2] += body.velocity[2] * dt;
		}

		this.state.timestamp += dt * 1000;
	}

	public setTimeScale(scale: number): void {
		this.timeScalar = scale;
	}
}
