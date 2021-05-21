import React from 'react';
import {mount, ReactWrapper} from 'enzyme';
import {stubUrl} from 'etna-js/spec/helpers';
import FlagsView from '../flags-view';

describe('FlagsView', () => {
  beforeEach(() => {
    const mockUsersResponse = [
      {
        email: 'janus@two-faces.org',
        name: 'Janus Bifrons',
        flags: ['inside', 'outside']
      },
      {
        email: 'portunus@two-faces.org',
        name: 'Portunus',
        flags: null
      }
    ];
    const mockProjectsResponse = [
      {project_name: 'doors', project_name_full: 'Oak doors'},
      {project_name: 'portals', project_name_full: 'Magic portals'}
    ];
    stubUrl({
      verb: 'get',
      url: 'https://janus.test/users',
      response: mockUsersResponse
    });
    stubUrl({
      verb: 'get',
      url: 'https://janus.test/projects',
      response: mockProjectsResponse
    });
  });

  it('renders correctly', () => {
    const component = mount(<FlagsView />);

    expect(component).toMatchSnapshot();
  });

  it('shows Add/Remove card when users are selected', () => {});

  it('filters users based on text search', () => {});

  it('filters users based on project search', () => {});

  it('filters users based on flag search', () => {});

  it('filters users based on all 3 search params', () => {});
});
