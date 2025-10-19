export enum MaterialType {
    Lambertian = 0,
    Metal = 1,
    Dielectric = 2,
}

export interface Material {
    albedo: [number, number, number];
    mat_type: MaterialType;
    fuzziness: number;
    refraction_index: number;
}

export interface Sphere {
    center: [number, number, number];
    radius: number;
    material: Material;
}

function rand(): number {
    return Math.random();
}

function randMinMax(min: number, max: number): number {
    return min + (max - min) * rand();
}

function randVec3(): [number, number, number] {
    return [rand(), rand(), rand()];
}

function randVec3MinMax(min: number, max: number): [number, number, number] {
    return [randMinMax(min, max), randMinMax(min, max), randMinMax(min, max)];
}

export function createRandomScene(): Sphere[] {
    const spheres: Sphere[] = [];

    // Ground
    spheres.push({
        center: [0.0, -1000.0, 0.0],
        radius: 1000.0,
        material: { albedo: [0.5, 0.5, 0.5], fuzziness: 0.0, refraction_index: 0.0, mat_type: MaterialType.Lambertian }
    });

    // Three large spheres
    spheres.push({ center: [0.0, 1.0, 0.0], radius: 1.0, material: { albedo: [1.0, 1.0, 1.0], fuzziness: 0.0, refraction_index: 1.5, mat_type: MaterialType.Dielectric } });
    spheres.push({ center: [-4.0, 1.0, 0.0], radius: 1.0, material: { albedo: [0.4, 0.2, 0.1], fuzziness: 0.0, refraction_index: 0.0, mat_type: MaterialType.Lambertian } });
    spheres.push({ center: [4.0, 1.0, 0.0], radius: 1.0, material: { albedo: [0.7, 0.6, 0.5], fuzziness: 0.0, refraction_index: 0.0, mat_type: MaterialType.Metal } });

    // Random small spheres up to 30 total
    for (let i = spheres.length; i < 30; i++) {
        const choose_mat = rand();
        const center: [number, number, number] = [
            randMinMax(-4.0, 4.0),
            0.2,
            randMinMax(-4.0, 4.0)
        ];

        const dist = Math.hypot(center[0] - 4.0, center[1] - 0.2, center[2] - 0.0);
        if (dist > 0.9) {
            if (choose_mat < 0.8) {
                const a = randVec3();
                const b = randVec3();
                const albedo: [number, number, number] = [a[0] * b[0], a[1] * b[1], a[2] * b[2]];
                spheres.push({ center, radius: 0.2, material: { albedo, fuzziness: 0.0, refraction_index: 0.0, mat_type: MaterialType.Lambertian } });
            } else if (choose_mat < 0.95) {
                const albedo = randVec3MinMax(0.5, 1.0);
                const fuzz = randMinMax(0.0, 0.5);
                spheres.push({ center, radius: 0.2, material: { albedo, fuzziness: fuzz, refraction_index: 0.0, mat_type: MaterialType.Metal } });
            } else {
                spheres.push({ center, radius: 0.2, material: { albedo: [1.0, 1.0, 1.0], fuzziness: 0.0, refraction_index: 1.5, mat_type: MaterialType.Dielectric } });
            }
        } else {
            spheres.push({ center: [0.0, 0.0, 0.0], radius: 0.1, material: { albedo: [1.0, 1.0, 1.0], fuzziness: 0.0, refraction_index: 0.0, mat_type: MaterialType.Lambertian } });
        }
    }

    return spheres;
}

export function serializeScene(spheres: Sphere[]): Float32Array {
    const floatsPerSphere = 12;
    const data = new Float32Array(spheres.length * floatsPerSphere);
    for (let i = 0; i < spheres.length; i++) {
        const s = spheres[i];
        const o = i * floatsPerSphere;
        data[o + 0] = s.center[0];
        data[o + 1] = s.center[1];
        data[o + 2] = s.center[2];
        data[o + 3] = s.radius;
        data[o + 4] = s.material.albedo[0];
        data[o + 5] = s.material.albedo[1];
        data[o + 6] = s.material.albedo[2];
        data[o + 7] = s.material.mat_type;
        data[o + 8] = s.material.fuzziness;
        data[o + 9] = s.material.refraction_index;
        data[o + 10] = 0.0; // padding
        data[o + 11] = 0.0; // padding
    }
    return data;
}


