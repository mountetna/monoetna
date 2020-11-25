import {projectNameFull} from '../janus';

describe('Janus Utils', () => {
  it('returns the matching project', () => {
    const state = {
      janus: {
        projects: [
          {
            project_name: 'xyz1',
            project_name_full: 'ACME corporation'
          }
        ]
      }
    };
    const result = projectNameFull(state, 'xyz1');

    expect(result).toEqual('ACME corporation');
  });

  it('returns the shortname if no matching project', () => {
    const state = {
      janus: {
        projects: [
          {
            project_name: 'xyz1',
            project_name_full: 'ACME corporation'
          }
        ]
      }
    };
    const result = projectNameFull(state, 'xyz2');

    expect(result).toEqual('xyz2');
  });
});
