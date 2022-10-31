import React from 'react';
import {Provider} from 'react-redux';
import {mockStore} from '../helpers';
import {render, screen, waitFor, fireEvent} from '@testing-library/react';
import '@testing-library/jest-dom/extend-expect';
import renderer from 'react-test-renderer';
import AttributeReport from '../../../lib/client/jsx/components/model_map/attribute_report';

const monster = require('../fixtures/template_monster.json');

describe('AttributeReport', () => {
  let store = mockStore({});

  beforeEach(() => {
    stubUrl({
      verb: 'post',
      path: '/retrieve',
      host: 'https://magma.test',
      status: 200,
      response: {}
    });

    stubUrl({
      verb: 'get',
      path: '/api/user/projects',
      host: 'https://janus.test',
      status: 200,
      response: {}
    });

    global.CONFIG = {
      magma_host: 'https://magma.test',
      janus_host: 'https://janus.test',
      project_name: 'labors'
    };
  });

  it('renders', () => {
    const tree = renderer
      .create(
        <Provider store={store}>
          <AttributeReport attribute={monster.attributes.name} />
        </Provider>
      )
      .toJSON();

    expect(tree).toMatchSnapshot();
  });

  it('renders with model action buttons for admin user', async () => {
    store = mockStore({
      magma: {
        models: {
          monster: require('../fixtures/template_monster.json')
        }
      },
      janus: {projects: require('../fixtures/project_names.json')},
      user: {
        permissions: {
          labors: {
            role: 'administrator'
          }
        }
      }
    });

    const {asFragment} = render(
      <Provider store={store}>
        <AttributeReport attribute={monster.attributes.name} />
      </Provider>
    );

    await waitFor(() => screen.getByText('Edit'));

    expect(screen.getByText('Edit')).toBeTruthy();
    expect(screen.queryByText('Remove')).toBeFalsy();
    expect(asFragment()).toMatchSnapshot();
  });

  it('renders with remove attribute buttons for admin user, removable attribute', async () => {
    store = mockStore({
      magma: {
        models: {
          monster: require('../fixtures/template_monster.json')
        }
      },
      janus: {projects: require('../fixtures/project_names.json')},
      user: {
        permissions: {
          labors: {
            role: 'administrator'
          }
        }
      }
    });

    const {asFragment} = render(
      <Provider store={store}>
        <AttributeReport attribute={monster.attributes.species} />
      </Provider>
    );

    await waitFor(() => screen.getByText('Remove'));

    expect(screen.getByText('Edit')).toBeTruthy();
    expect(screen.getByText('Remove')).toBeTruthy();
    expect(asFragment()).toMatchSnapshot();
  });
});
