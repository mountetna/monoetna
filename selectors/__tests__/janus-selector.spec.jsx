import {selectProjects} from '../janus-selector';

describe('Janus Selector', () => {
  it('returns projects', () => {
    const projects = [{project_name: 'xyz1'}];
    const result = selectProjects({
      janus: {projects}
    });

    expect(result).toEqual(projects);
  });
});
