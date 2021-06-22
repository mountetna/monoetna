import {projectNameFull, parseToken} from '../janus';

describe('Janus Utils', () => {
  describe('projectNameFull', () => {
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

  describe('parseToken', () => {
    function generateTestToken(params) {
      return `header.${btoa(JSON.stringify(params))}.signature`;
    }

    it('correctly parses token for user with flags', () => {
      const userInfo = {
        email: 'janus@two-faces.org',
        name: 'Janus Portunus',
        perm: 'a:project1,project2;v:project3',
        flags: 'red;yellow;blue'
      };
      const result = parseToken(generateTestToken(userInfo));

      expect(result).toEqual({
        email: 'janus@two-faces.org',
        name: 'Janus Portunus',
        permissions: {
          project1: {
            privileged: false,
            project_name: 'project1',
            role: 'administrator'
          },
          project2: {
            privileged: false,
            project_name: 'project2',
            role: 'administrator'
          },
          project3: {
            privileged: false,
            project_name: 'project3',
            role: 'viewer'
          }
        },
        flags: ['red', 'yellow', 'blue']
      });
    });

    it('correctly parses token for user with no flags', () => {
      const userInfo = {
        email: 'janus@two-faces.org',
        name: 'Janus Portunus',
        perm: 'a:project1,project2;v:project3'
      };

      const result = parseToken(generateTestToken(userInfo));

      expect(result).toEqual({
        email: 'janus@two-faces.org',
        name: 'Janus Portunus',
        permissions: {
          project1: {
            privileged: false,
            project_name: 'project1',
            role: 'administrator'
          },
          project2: {
            privileged: false,
            project_name: 'project2',
            role: 'administrator'
          },
          project3: {
            privileged: false,
            project_name: 'project3',
            role: 'viewer'
          }
        },
        flags: []
      });
    });
  });
});
