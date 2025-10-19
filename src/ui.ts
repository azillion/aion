import { GUI } from 'dat.gui';
import type { Renderer } from './renderer';
import { spectralResponses } from './spectral';
import { themes } from './theme';

export class UI {
    private renderer: Renderer;
    private gui: GUI;
    private settings = {
        'Time Scale (days/sec)': 30,
        'Focus': 'Sun',
        'Theme': 'amber',
        'Sensor': 'Visible (Y)',
    };

    constructor(renderer: Renderer) {
        this.renderer = renderer;
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


