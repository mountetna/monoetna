import React from 'react';
import {Provider} from 'react-redux';
import {mount, shallow} from 'enzyme';
import renderer from 'react-test-renderer';
import {mockStore} from '../../spec/helpers';
import RootView from '../RootView';

const permissions = require('../../spec/fixtures/home_page_permissions.json');
const projects = require('../../spec/fixtures/project_names.json')

describe('RootView', () => {
  let store;

  beforeEach(() => {
    store = mockStore({
      user: { permissions },
      janus: { projects }
    });
  });
  it('renders', () => {
    const tree = renderer
    .create(<Provider store={ store }>
      <RootView />
    </Provider>)
      .toJSON();

    expect(tree).toMatchSnapshot();
  });
});
