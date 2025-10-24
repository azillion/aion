import { GUI } from 'dat.gui';
import type { Renderer } from './renderer/index';
import { spectralResponses } from './spectral';
import { themes } from './theme';
import { AppState, CameraMode, ReferenceFrame } from './state';
import { G } from './shared/constants';
import type { Authority } from './authority/authority';
import { CameraManager } from './camera/manager';

export class UI {
    private renderer: Renderer;
    private state: AppState;
    private authority: Authority;
    private gui: GUI;
    private focusController: dat.GUIController | null = null;
    private cameraManager: CameraManager;
    private timeFolder: GUI;
    private settings = {
        'Focus': 'Sun',
        'Theme': 'white',
        'Sensor': 'Full Color',
        'Show Orbits': false,
        'Add Asteroid': () => this.addAsteroid(),
        'Warp x1': () => this.setWarp(1),
        'Warp x100': () => this.setWarp(100),
        'Warp x1k': () => this.setWarp(1000),
        'Warp x10k': () => this.setWarp(10000),
        'Warp x100k': () => this.setWarp(100000),
        'Warp x1M': () => this.setWarp(1000000),
        'Time Scale': 1.0, // We keep the slider for variable control
        'Current Warp': 'x1',
        'Toggle View': () => {
            if (this.state.cameraMode === CameraMode.SYSTEM_MAP) {
                this.state.cameraMode = CameraMode.SHIP_RELATIVE;
            } else if (this.state.cameraMode === CameraMode.SHIP_RELATIVE) {
                this.state.cameraMode = CameraMode.SYSTEM_MAP;
            } else {
                this.state.cameraMode = CameraMode.SHIP_RELATIVE;
            }
        }
    };

    constructor(renderer: Renderer, state: AppState, authority: Authority, cameraManager: CameraManager) {
        this.renderer = renderer;
        this.state = state;
        this.authority = authority;
        this.cameraManager = cameraManager;
        this.gui = new GUI();

        this.timeFolder = this.gui.addFolder('Time Control');
        this.timeFolder.add(this.settings, 'Warp x1');
        this.timeFolder.add(this.settings, 'Warp x100');
        this.timeFolder.add(this.settings, 'Warp x1k');
        this.timeFolder.add(this.settings, 'Warp x10k');
        this.timeFolder.add(this.settings, 'Warp x100k');
        this.timeFolder.add(this.settings, 'Warp x1M');
        this.timeFolder.add(this.settings, 'Current Warp').listen();
        this.timeFolder
            .add(this.settings, 'Time Scale', 0, 10000)
            .onChange((value: number) => this.setWarp(value))
            .listen();
        this.timeFolder.open();

        this._createSceneControls();

        this.gui.add(this.settings, 'Add Asteroid');

        this.gui.add(this.settings, 'Theme', Object.keys(themes))
            .onChange(() => this.updateTheme());

        this.gui.add(this.settings, 'Sensor', Object.keys(spectralResponses))
            .onChange(() => this.updateTheme());

        this.gui.add(this.settings, 'Show Orbits')
            .onChange((value: boolean) => {
                this.state.showOrbits = value;
                if (!value) {
                    this.renderer.clearOrbitHistory();
                }
            });

        this.gui.add(this.state, 'crtIntensity', 0.0, 1.0).name('CRT Effect');

        this.gui.add(this.state, 'debugTierView', {
            'Normal': -1,
            'Near Tier Only': 0,
            'Mid Tier Only': 1,
            'Far Tier Only': 2,
        }).name('Debug View');

        this.gui.add(this.settings, 'Toggle View');

        this.gui.add(this.state, 'referenceFrame', Object.values(ReferenceFrame)).name('Reference Frame');
    }

    public setFocus(bodyName: string) {
        this.settings['Focus'] = bodyName;
        if (this.focusController) {
            this.focusController.updateDisplay();
        }
    }

    private setWarp(scale: number) {
        this.authority.setTimeScale(scale);
        this.settings['Time Scale'] = scale;
        const scaleNames: {[key: number]: string} = {1: 'x1', 100: 'x100', 1000: 'x1k', 10000: 'x10k', 100000: 'x100k', 1000000: 'x1M'};
        this.settings['Current Warp'] = scaleNames[scale] || `${scale.toPrecision(2)}x`;
    }

    private updateTheme = () => {
        this.renderer.setTheme(this.settings['Theme'], this.settings['Sensor']);
    }

    private async _createSceneControls() {
        if (this.focusController) {
            this.gui.remove(this.focusController);
        }
        const initialState = await this.authority.query();
        const bodyNames = initialState.bodies.map(b => b.name);
        const bodyIds = initialState.bodies.map(b => b.id);

        const playerIndex = bodyIds.indexOf('player-ship');
        const cam = this.cameraManager.getCamera();
        if (playerIndex !== -1) {
            this.settings['Focus'] = 'AION-1';
            cam.focusBodyId = 'player-ship';
            cam.pendingFrame = true;
        } else {
            this.settings['Focus'] = bodyNames[0];
            cam.focusBodyId = bodyIds[0];
            cam.pendingFrame = true;
        }

        this.focusController = this.gui.add(this.settings, 'Focus', bodyNames)
            .onChange((selectedName: string) => {
                const selectedIndex = bodyNames.indexOf(selectedName);
                const selectedId = bodyIds[selectedIndex];
                const cam2 = this.cameraManager.getCamera();
                cam2.focusBodyId = selectedId;
                cam2.pendingFrame = true;
            });
    }

    private async addAsteroid() {
        const systemState = await this.authority.query();
        const bodies = systemState.bodies;
        if (bodies.length === 0) return;
        const focusId = this.cameraManager.getCamera().focusBodyId;
        const parent = (focusId ? bodies.find(b => b.id === focusId) : null) || bodies[0];

        const orbitDistance = Math.max(1, parent.radius * 50);
        const position = [
            parent.position[0] + orbitDistance,
            parent.position[1],
            parent.position[2],
        ] as [number, number, number];
        const orbitalSpeed = Math.sqrt(Math.max(0, G * parent.mass / orbitDistance));
        const velocity = [
            parent.velocity[0],
            parent.velocity[1] + orbitalSpeed,
            parent.velocity[2],
        ] as [number, number, number];

        const name = `Asteroid-${Date.now().toString(36)}`;
        const radius = 10;
        const mass = 1e15;
        const albedo: [number, number, number] = [0.6, 0.6, 0.6];

        this.authority.addBody({ name, position, velocity, radius, mass, albedo });

        setTimeout(() => this._createSceneControls(), 100);
    }
}
