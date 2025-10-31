import { WasmAuthority } from './wasmAuthority';
import type { ClientToServerMessage, ServerToClientMessage } from '@shared/messages';
import type { Body, Ship, SystemState } from '@shared/types';
import { G } from '@shared/constants';

const authority = new WasmAuthority();

async function init(): Promise<void> {
  // Build initial system state (ported from LocalAuthority constructor)
  const sun: Body = {
    id: 'sol', name: 'Sun', position: [0, 0, 0], velocity: [0, 0, 0], radius: 696340, mass: 1.989e30, albedo: [1.0, 1.0, 0.9], emissive: [1000, 1000, 1000]
  } as any;
  const earth_dist = 149.6e6;
  const earth_vel = Math.sqrt(G * sun.mass / earth_dist);
  const earth: Body = {
    id: 'earth', name: 'Earth', position: [0, earth_dist, 0], velocity: [-earth_vel, 0, 0], radius: 6371, mass: 5.972e24, albedo: [0.2, 0.3, 0.8]
  } as any;
  const moon_dist = 384400;
  const moon_vel = Math.sqrt(G * earth.mass / moon_dist);
  const moon: Body = {
    id: 'moon', name: 'Moon', position: [0, earth_dist - moon_dist, 0], velocity: [-earth_vel + moon_vel, 0, 0], radius: 1737, mass: 7.347e22, albedo: [0.5, 0.5, 0.5]
  } as any;
  const playerShip: Ship = {
    id: 'player-ship', name: 'AION-1',
    position: [earth.position[0] + (earth.radius + 35786), earth.position[1], earth.position[2]],
    velocity: [earth.velocity[0], earth.velocity[1] + Math.sqrt(G * earth.mass / (earth.radius + 35786)), earth.velocity[2]],
    radius: 0.1, mass: 1e6, albedo: [0.8, 0.8, 0.9], emissive: [0.1, 0.3, 1.0],
    orientation: [0, 0, 0, 1], angularVelocity: [0, 0, 0], thrust: [0, 0, 0]
  } as any;

  const initial: SystemState = { timestamp: Date.now(), bodies: [sun, earth, moon, playerShip] };
  await authority.initialize(initial);
}

onmessage = async (e: MessageEvent<ClientToServerMessage>) => {
  const msg = e.data;
  switch (msg.type) {
    case 'query': {
      const state = await authority.query();
      const reply: ServerToClientMessage = { type: 'queryResult', queryId: msg.queryId, state };
      postMessage(reply);
      break;
    }
    case 'tick': {
      await authority.tick(msg.deltaTime, msg.input as any);
      break;
    }
    case 'setTimeScale': {
      await authority.setTimeScale(msg.scale);
      break;
    }
    case 'addBody': {
      await authority.addBody(msg.body);
      break;
    }
    case 'autoLand': {
      await authority.autoLand(msg.targetBodyId);
      break;
    }
    case 'teleportToSurface': {
      await authority.teleportToSurface(msg.targetBodyId);
      break;
    }
  }
};

void init();


