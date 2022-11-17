import React from 'react';
import {Provider} from 'react-redux';
import {mockStore, stubUrl} from '../helpers';
import {render, screen, waitFor, fireEvent} from '@testing-library/react';
import '@testing-library/jest-dom/extend-expect';
import renderer from 'react-test-renderer';
import AttributeReport from '../../../lib/client/jsx/components/model_map/attribute_report';

const monster = require('../fixtures/template_monster.json');
const habitat = require('../fixtures/template_habitat.json');

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

  it('renders without model action buttons for editor user', async () => {
    store = mockStore({
      magma: {
        models: {
          monster: {
            template: require('../fixtures/template_monster.json')
          },
          habitat: {
            template: require('../fixtures/template_habitat.json')
          }
        }
      },
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
        <AttributeReport attribute={monster.attributes.species} />
      </Provider>
    );

    await waitFor(() => screen.getByText('Species'));

    expect(screen.queryByText('Edit')).toBeFalsy();
    expect(screen.queryByText('Remove')).toBeFalsy();
    expect(asFragment()).toMatchSnapshot();
  });

  it('renders without model action buttons for admin user, non-editable attribute', async () => {
    store = mockStore({
      magma: {
        models: {
          monster: {
            template: require('../fixtures/template_monster.json')
          },
          habitat: {
            template: require('../fixtures/template_habitat.json')
          }
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

    await waitFor(() => screen.getByText('Name'));

    expect(screen.queryByText('Edit')).toBeFalsy();
    expect(screen.queryByText('Remove')).toBeFalsy();
    expect(asFragment()).toMatchSnapshot();
  });

  it('renders with remove attribute buttons for admin user, removable attribute', async () => {
    store = mockStore({
      magma: {
        models: {
          monster: {
            template: require('../fixtures/template_monster.json')
          },
          habitat: {
            template: require('../fixtures/template_habitat.json')
          }
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
        <AttributeReport
          attribute={monster.attributes.species}
          isAdminUser={true}
        />
      </Provider>
    );

    await waitFor(() => screen.getByText('Remove'));

    expect(screen.getByText('Edit')).toBeTruthy();
    expect(screen.getByText('Remove')).toBeTruthy();
    expect(screen.queryByText('Remove Link')).toBeFalsy();
    expect(asFragment()).toMatchSnapshot();
  });

  it('renders with remove link button for admin user, link attribute', async () => {
    store = mockStore({
      magma: {
        models: {
          monster: {
            template: require('../fixtures/template_monster.json')
          },
          habitat: {
            template: require('../fixtures/template_habitat.json')
          }
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
        <AttributeReport
          attribute={monster.attributes.habitat}
          isAdminUser={true}
        />
      </Provider>
    );

    await waitFor(() => screen.getByText('Remove Link'));

    expect(screen.queryByText('Edit')).toBeFalsy();
    expect(screen.getByText('Remove Link')).toBeTruthy();
    expect(asFragment()).toMatchSnapshot();
  });

  it('renders with remove link button for admin user, collection (link) attribute', async () => {
    store = mockStore({
      magma: {
        models: {
          monster: {
            template: require('../fixtures/template_monster.json')
          },
          habitat: {
            template: require('../fixtures/template_habitat.json')
          }
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
        <AttributeReport
          attribute={habitat.attributes.monster}
          isAdminUser={true}
        />
      </Provider>
    );

    await waitFor(() => screen.getByText('Remove Link'));

    expect(screen.queryByText('Edit')).toBeFalsy();
    expect(screen.getByText('Remove Link')).toBeTruthy();
    expect(asFragment()).toMatchSnapshot();
  });
});
