import React from 'react';
import {Provider} from 'react-redux';
import { act, create } from 'react-test-renderer';
import {mockStore} from '../test-helpers';
import {stubUrl} from 'etna-js/spec/helpers';
import ProjectView from '../project-view';

describe('ProjectView', () => {
  const setRole = (role) => {
    stubUrl({
      verb: 'get',
      path: '/api/admin/ports/info',
      status: 200,
      response: {
        project: {
          permissions: [
            {
              affiliation: null,
              privileged: true,
              project_name: "ports",
              role,
              user_email: "janus@two-faces.org",
              user_name: "Janus Bifrons"
            },
            {
              affiliation: "ILWU Local 34",
              privileged: null,
              project_name: "ports",
              role: "editor",
              user_email: "portunus@two-faces.org",
              user_name: "Portunus"
            }
          ],
          project_name: "ports",
          project_name_full: "Ports"
        }
      }
    });

    return mockStore({
      user: {
        email: "janus@two-faces.org",
        name: "Janus Bifrons",
        permissions: {
          ports: {
            privileged: true,
            project_name: 'ports',
            role
          }
        }
      }
    });
  }

  it('renders statically for an editor', async () => {
    let store = setRole('editor');

    let component = create(
      <Provider store={store}>
        <ProjectView project_name='ports'/>
      </Provider>
    );

    await act( async () => {
        await new Promise((resolve) => setTimeout(resolve, 20));
    });

    expect( component.toJSON() ).toMatchSnapshot()
  });

  it('renders with inputs for an administrator', async () => {
    let store = setRole('administrator');

    let component = create(
      <Provider store={store}>
        <ProjectView project_name='ports'/>
      </Provider>
    );

    await act( async () => {
        await new Promise((resolve) => setTimeout(resolve, 20));
    });

    expect( component.toJSON() ).toMatchSnapshot()
  });
});
