import React from 'react';
import {Provider} from 'react-redux';
import {render, screen, waitFor, fireEvent} from '@testing-library/react';
import '@testing-library/jest-dom/extend-expect';
import renderer from 'react-test-renderer';
import {mockStore, stubUrl} from '../helpers';
import ModelMap from '../../../lib/client/jsx/components/model_map';

const models = {
  monster: {template: require('../fixtures/template_monster.json')},
  labor: {template: require('../fixtures/template_labor.json')},
  project: {template: require('../fixtures/template_project.json')}
};

describe('ModelMap for model with no records', () => {
  let store;

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

    stubUrl({
      verb: 'post',
      path: '/query',
      host: 'https://magma.test',
      status: 200,
      response: {
        answer: 0
      },
      request: () => true,
      times: 20
    });

    global.CONFIG = {
      magma_host: 'https://magma.test',
      janus_host: 'https://janus.test',
      project_name: 'labors'
    };
  });

  it('renders', () => {
    store = mockStore({
      magma: {models},
      janus: {projects: require('../fixtures/project_names.json')}
    });

    // Wrap with Provider here so store gets passed down to child components in Context
    const tree = renderer
      .create(
        <Provider store={store}>
          <ModelMap />
        </Provider>
      )
      .toJSON();

    expect(tree).toMatchSnapshot();
  });

  it('renders with model action buttons for admin user (disabled reparent model btn while awaiting record count)', async () => {

    store = mockStore({
      magma: {models},
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
        <ModelMap />
      </Provider>
    );

    await waitFor(() => screen.findByText('Attribute'));

    expect(screen.getByTitle('Add Link')).toBeTruthy();
    expect(screen.getByTitle('Add Attribute')).toBeTruthy();
    expect(screen.getByTitle('Add Model')).toBeTruthy();
    expect(screen.getByTitle('Determining if reparenting is possible')).toBeTruthy();
    expect(screen.getByTitle('Cannot remove a model with any attribute of the types: table, child, collection')).toBeTruthy();
    expect(asFragment()).toMatchSnapshot();
  });

  it('renders with action buttons and reparent model enabled for admin user', async () => {  
    store = mockStore({
      magma: {models},
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
        <ModelMap />
      </Provider>
    );

    await waitFor(() => screen.findByText('Attribute'));

    const monsterModel = screen.getByText('monster');
    fireEvent.click(monsterModel);

    await waitFor(() => screen.getByText('0 records'));
    expect(screen.getByTitle('Add Link')).toBeTruthy();
    expect(screen.getByTitle('Add Attribute')).toBeTruthy();
    expect(screen.getByTitle('Add Model')).toBeTruthy();
    expect(screen.getByTitle('Reparent Model')).toBeTruthy();
    expect(asFragment()).toMatchSnapshot();
  });

  it('renders with remove model enabled for leaf model and admin user', async () => {
    store = mockStore({
      magma: {models: {...models, aspect: {template: require('../fixtures/template_aspect.json')}}},
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
        <ModelMap />
      </Provider>
    );

    await waitFor(() => screen.findByText('Attribute'));

    const aspectModel = screen.getByText('aspect');
    fireEvent.click(aspectModel);

    expect(screen.getByTitle('Add Link')).toBeTruthy();
    expect(screen.getByTitle('Add Attribute')).toBeTruthy();
    expect(screen.getByTitle('Add Model')).toBeTruthy();
    expect(screen.getByTitle('Remove Model')).toBeTruthy();
    expect(asFragment()).toMatchSnapshot();
  });

  it('renders without model action buttons for editor user', async () => {
    store = mockStore({
      magma: {models},
      janus: {projects: require('../fixtures/project_names.json')},
      user: {
        permissions: {
          labors: {
            role: 'editor'
          }
        }
      }
    });

    const {asFragment} = render(
      <Provider store={store}>
        <ModelMap />
      </Provider>
    );

    await waitFor(() => screen.findByText('monster'));

    expect(screen.queryByTitle('Add Link')).toBeFalsy();
    expect(screen.queryByTitle('Add Attribute')).toBeFalsy();
    expect(screen.queryByTitle('Add Model')).toBeFalsy();
    expect(screen.queryByText('Reparent Model')).toBeFalsy();
    expect(asFragment()).toMatchSnapshot();
  });

  it('renders an incomplete graph', () => {
    let {monster, project} = models;

    store = mockStore({
      magma: {models: {monster, project}},
      janus: {projects: []}
    });

    // Wrap with Provider here so store gets passed down to child components in Context
    const tree = renderer
      .create(
        <Provider store={store}>
          <ModelMap />
        </Provider>
      )
      .toJSON();

    expect(tree).toMatchSnapshot();
  });
});
