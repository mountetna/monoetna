import React from 'react';
import {Provider} from 'react-redux';
import { act, create } from 'react-test-renderer';
import {mockStore} from '../test-helpers';
import {stubUrl} from 'etna-js/spec/helpers';
import JanusAdmin from '../janus-admin';

describe('JanusAdmin', () => {
  let store;

  beforeEach(() => {
    store = mockStore({
      user: {
	email: "janus@two-faces.org",
	name: "Janus Bifrons",
        permissions: { 
          administration: {
            privileged: false,
            project_name: 'administration',
            role: 'editor'
          }
        }
      }
    });
  });

  it('renders', async () => {
    const initialStubs = [
      stubUrl({
        verb: 'get',
        path: '/allprojects',
        status: 200,
        response: {
          projects: [
            { project_name: "gateway", project_name_full: "Gateway"},
            { project_name: "tunnel", project_name_full: "Tunnel"},
            { project_name: "mirror", project_name_full: "Mirror"}
          ]
        }
      })
    ];

    let component = create(
      <Provider store={store}>
        <JanusAdmin/>
      </Provider>
    );

    await act( async () => {
        await new Promise((resolve) => setTimeout(resolve, 15));
    });

    expect( component.toJSON() ).toMatchSnapshot()
  });
});
