import { GUI } from 'dat.gui';
import type { Authority } from './authority/authority';
import type { Camera } from './camera';

export class UI {
	private authority: Authority;
	private camera: Camera;
	private gui: GUI;
	private settings = {
		'Time Scale (days/sec)': 30,
		'Focus': 'Sun',
	};

	constructor(authority: Authority, camera: Camera) {
		this.authority = authority;
		this.camera = camera;
		this.gui = new GUI();

		this.gui
			.add(this.settings, 'Time Scale (days/sec)', 0, 365)
			.onChange((value: number) => {
				const secondsPerDay = 24 * 60 * 60;
				this.authority.setTimeScale(value * secondsPerDay);
			});


		this._createSceneControls();
	}

	private async _createSceneControls() {
		const initialState = await this.authority.query();
		const bodyNames = initialState.bodies.map(b => b.name);
		const bodyIds = initialState.bodies.map(b => b.id);

		this.settings['Focus'] = bodyNames[0];
		this.camera.setFocus(bodyIds[0]);

		this.gui.add(this.settings, 'Focus', bodyNames)
			.onChange((selectedName: string) => {
				const selectedIndex = bodyNames.indexOf(selectedName);
				const selectedId = bodyIds[selectedIndex];
				this.camera.setFocus(selectedId);
			});
	}
}


