import { GUI } from 'dat.gui';
import type { Renderer } from './renderer';
import { spectralResponses } from './spectral';
import { themes } from './theme';
import { AppState, ViewMode } from './state';

export class UI {
    private renderer: Renderer;
    private state: AppState;
    private gui: GUI;
    private settings = {
        'Time Scale (days/sec)': 30,
        'Focus': 'Sun',
        'Theme': 'amber',
        'Sensor': 'Visible (Y)',
        'Toggle View': () => {
            if (this.state.viewMode === ViewMode.System) {
                this.state.viewMode = ViewMode.Galaxy;
                this.renderer.camera.eye = [0, 0, 150];
                this.renderer.camera.look_at = [0, 0, 0];
            } else {
                this.state.viewMode = ViewMode.System;
                this.renderer.camera.eye = [0, 0.3, 1.0];
                this.renderer.camera.look_at = [0, 0, -2.0];
            }
        }
    };

    constructor(renderer: Renderer, state: AppState) {
        this.renderer = renderer;
        this.state = state;
        this.gui = new GUI();

        this.gui
            .add(this.settings, 'Time Scale (days/sec)', 0, 365)
            .onChange((value: number) => {
                const secondsPerDay = 24 * 60 * 60;
                this.renderer.authority.setTimeScale(value * secondsPerDay);
            });

        this._createSceneControls();

        this.gui.add(this.settings, 'Theme', Object.keys(themes))
            .onChange(() => this.updateTheme());

        this.gui.add(this.settings, 'Sensor', Object.keys(spectralResponses))
            .onChange(() => this.updateTheme());

        this.gui.add(this.settings, 'Toggle View');
    }

    private updateTheme = () => {
        this.renderer.setTheme(this.settings['Theme'], this.settings['Sensor']);
    }

    private async _createSceneControls() {
        const initialState = await this.renderer.authority.query();
        const bodyNames = initialState.bodies.map(b => b.name);
        const bodyIds = initialState.bodies.map(b => b.id);

        this.settings['Focus'] = bodyNames[0];
        this.renderer.camera.setFocus(bodyIds[0]);

        this.gui.add(this.settings, 'Focus', bodyNames)
            .onChange((selectedName: string) => {
                const selectedIndex = bodyNames.indexOf(selectedName);
                const selectedId = bodyIds[selectedIndex];
                this.renderer.camera.setFocus(selectedId);
            });
    }
}


