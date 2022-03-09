import React from 'react';
import {Provider} from 'react-redux';
import { act, create } from 'react-test-renderer';
import {mockStore} from '../test-helpers';
import {stubUrl} from 'etna-js/spec/helpers';
import JanusMain from '../janus-main';

describe('JanusMain', () => {
  let store;

  beforeEach(() => {
    store = mockStore({
      user: {
	email: "janus@two-faces.org",
	name: "Janus Bifrons",
	permissions: { }
      }
    });
  });

  it('renders', async () => {
    const initialStubs = [
      stubUrl({
        verb: 'get',
	path: '/api/user/projects',
        status: 200,
        response: {projects: [
          {
            project_name: "tunnel",
            project_name_full: "Tunnel",
            role: "viewer",
            privileged: true
          }, {
            project_name: "mirror",
            project_name_full: "Mirror",
            role: "editor",
            privileged: null
          }, {
            project_name: "gateway",
            project_name_full: "Gateway",
            role: "editor",
            privileged: null
          } ]}
      })
    ];

    let component = create(
      <Provider store={store}>
        <JanusMain/>
      </Provider>
    );

    await act( async () => {
        await new Promise((resolve) => setTimeout(resolve, 15));
    });

    expect( component.toJSON() ).toMatchSnapshot()
  });
});
