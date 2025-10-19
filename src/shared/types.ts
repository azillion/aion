export type Vec3 = [number, number, number];

export interface Body {
	id: string;
	name: string;
	position: Vec3;
	velocity: Vec3;
	radius: number;
	mass: number;
	albedo: Vec3;
	emissive?: Vec3;
}

export interface SystemState {
	timestamp: number;
	bodies: Body[];
}
