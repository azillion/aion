import type { Authority } from './authority';
import type { Body, SystemState } from '../shared/types';

export class LocalAuthority implements Authority {
	private state: SystemState;

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

		const earth: Body = {
			id: 'earth',
			name: 'Earth',
			position: [149.6e6, 0, 0],
			velocity: [0, 29.78, 0],
			radius: 6371,
			mass: 5.972e24,
			albedo: [0.2, 0.3, 0.8],
		};
		
		const moon: Body = {
			id: 'moon',
			name: 'Moon',
			position: [149.6e6 + 384400, 0, 0],
			velocity: [0, 29.78 + 1.022, 0],
			radius: 1737,
			mass: 7.347e22,
			albedo: [0.5, 0.5, 0.5],
		};

		this.state = {
			timestamp: Date.now(),
			bodies: [sun, earth, moon],
		};
	}

	async query(): Promise<SystemState> {
		return this.state;
	}

	async tick(deltaTime: number): Promise<void> {
		const timeScalar = 1_000_000;
		const dt = deltaTime * timeScalar;

		this.state.bodies.forEach(body => {
			body.position[0] += body.velocity[0] * dt;
			body.position[1] += body.velocity[1] * dt;
			body.position[2] += body.velocity[2] * dt;
		});

		this.state.timestamp += dt * 1000;
	}
}
