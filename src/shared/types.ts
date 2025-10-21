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

export interface Ship extends Body {
	orientation: [number, number, number, number];
	thrust: Vec3;
}

export interface SystemState {
	timestamp: number;
	bodies: Body[];
}

export interface Theme {
  bg: Vec3;
  fg: Vec3;
  accent: Vec3;
}

export interface Star {
  position: Vec3;
  color: Vec3;
  size: number;
}
