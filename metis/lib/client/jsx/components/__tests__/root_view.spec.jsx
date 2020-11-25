import React from 'react';
import {Provider} from 'react-redux';
import renderer from 'react-test-renderer';
import thunk from 'redux-thunk';
import configureMockStore from 'redux-mock-store';
import RootView from '../root-view';

const permissions = require('./fixtures/root_view_permissions.json');
const projects = require('./fixtures/project_names.json');

describe('RootView', () => {
  it('renders', () => {
    const mockStore = configureMockStore([thunk]);

    const tree = renderer
      .create(
        <Provider
          store={mockStore({
            user: {
              permissions: permissions
            },
            janus: {
              projects: projects
            }
          })}
        >
          <RootView />
        </Provider>
      )
      .toJSON();

    expect(tree).toMatchSnapshot();
  });
});
