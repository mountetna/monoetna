import {projectNameFull} from '../janus';

describe('Janus Utils', () => {
  it('returns the matching project', () => {
    const projects = [
      {
        project_name: 'xyz1',
        project_name_full: 'ACME corporation'
      }
    ];
    const result = projectNameFull(projects, 'xyz1');

    expect(result).toEqual('ACME corporation');
  });

  it('returns null if no matching project', () => {
    const projects = [
      {
        project_name: 'xyz1',
        project_name_full: 'ACME corporation'
      }
    ];
    const result = projectNameFull(projects, 'xyz2');

    expect(result).toEqual(null);
  });

  it('returns null if no projects', () => {
    const result = projectNameFull(null, 'xyz2');

    expect(result).toEqual(null);
  });
});
