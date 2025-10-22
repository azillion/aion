import { GUI } from 'dat.gui';
import type { Renderer } from './renderer/index';
import { spectralResponses } from './spectral';
import { themes } from './theme';
import { AppState, CameraMode, ReferenceFrame } from './state';
import { G } from './shared/constants';

export class UI {
    private renderer: Renderer;
    private state: AppState;
    private gui: GUI;
    private focusController: dat.GUIController | null = null;
    private settings = {
        'Time Scale': 1.0,
        'Focus': 'Sun',
        'Theme': 'white',
        'Sensor': 'Visible (Y)',
        'Show Orbits': false,
        'Add Asteroid': () => this.addAsteroid(),
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

    constructor(renderer: Renderer, state: AppState) {
        this.renderer = renderer;
        this.state = state;
        this.gui = new GUI();

        this.gui
            .add(this.settings, 'Time Scale', 0, 10000)
            .onChange((value: number) => {
                this.renderer.authority.setTimeScale(value);
            });

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

        this.gui.add(this.settings, 'Toggle View');

        this.gui.add(this.state, 'referenceFrame', Object.values(ReferenceFrame)).name('Reference Frame');
    }

    public setFocus(bodyName: string) {
        this.settings['Focus'] = bodyName;
        if (this.focusController) {
            this.focusController.updateDisplay();
        }
    }

    private updateTheme = () => {
        this.renderer.setTheme(this.settings['Theme'], this.settings['Sensor']);
    }

    private async _createSceneControls() {
        if (this.focusController) {
            this.gui.remove(this.focusController);
        }
        const initialState = await this.renderer.authority.query();
        const bodyNames = initialState.bodies.map(b => b.name);
        const bodyIds = initialState.bodies.map(b => b.id);

        const playerIndex = bodyIds.indexOf('player-ship');
        const cam = this.renderer.getCamera();
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
                const cam2 = this.renderer.getCamera();
                cam2.focusBodyId = selectedId;
                cam2.pendingFrame = true;
            });
    }

    private async addAsteroid() {
        const systemState = await this.renderer.authority.query();
        const bodies = systemState.bodies;
        if (bodies.length === 0) return;
        const focusId = this.renderer.getCamera().focusBodyId;
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

        this.renderer.authority.addBody({ name, position, velocity, radius, mass, albedo });

        setTimeout(() => this._createSceneControls(), 100);
    }
}
