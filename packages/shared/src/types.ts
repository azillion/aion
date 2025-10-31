export type Vec3 = [number, number, number];

// Atmospheric gas composition by fraction (should sum to ~1.0)
export interface AtmosphereComposition {
    [gas: string]: number;
}

export interface TerrainParams {
    radius: number;
    seaLevel: number;
    maxHeight: number;
    noiseSeed: number;
    atmosphere?: AtmosphereComposition;
}

export interface Body {
	id: string;
	name: string;
	position: Vec3;
	velocity: Vec3;
	radius: number;
	mass: number;
	albedo: Vec3;
	emissive?: Vec3 | null;
    terrain?: TerrainParams;
}

export interface Ship extends Body {
	orientation: [number, number, number, number];
	angularVelocity: Vec3;
	thrust: Vec3;
}

export interface SystemState {
	timestamp: number;
	bodies: Body[];
	ships: Ship[];
	flags?: { precision?: boolean; killRotation?: boolean };
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

export interface Orbit {
	semiMajorAxis: number;
	eccentricity: number;
	// Basis vectors for the orbital plane
	p: Vec3;
	q: Vec3;
	// Standard gravitational parameter
	mu: number;
	// True anomaly at the time the orbit was calculated (epoch)
	trueAnomalyAtEpoch: number;
	// Simulation timestamp of the epoch (seconds)
	epoch: number;
}

export interface FrameData {
	rawState: SystemState;
	bodiesToRender: Body[];
	camera: any;
	systemScale: number;
	viewport: { width: number, height: number };
	deltaTime: number;
	cameraMode: number;
	playerShipId: string | null;
  	dominantLight?: Body;
  	worldCameraEye?: Vec3;
    debugTierView?: number;
    showOrbits: boolean;
    showAtmosphere: boolean;
    unscaledBodiesForMap?: Body[];
}
