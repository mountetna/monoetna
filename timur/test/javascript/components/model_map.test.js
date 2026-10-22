import React from 'react';
import {Provider} from 'react-redux';
import {render, screen, waitFor, fireEvent, within} from '@testing-library/react';
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

  it('runs and displays the project audit returned by the API', async () => {
    stubUrl({
      verb: 'get',
      path: '/cross_project_search_audit/labors',
      host: 'https://magma.test',
      status: 200,
      response: {
        project: 'labors',
        conforming: false,
        template_mapping: {
          mapped_models: [{model: 'labor', template_model: 'labor'}],
          unmapped_models: ['monster'],
          invalid_mappings: [
            {model: 'prize', template_model: 'award'},
            {model: 'cup', template_model: null}
          ],
          required_column_results: [
            {model: 'labor', template_model: 'labor', column: 'name', present: true},
            {model: 'labor', template_model: 'labor', column: 'species', present: false}
          ],
          ontology_value_results: [
            {model: 'labor', template_model: 'labor', column: 'organism', value: 'PBMC', valid: false}
          ],
          ontology_valid_count: 1,
          ontology_mismatches: [{
            model: 'labor',
            column: 'tissue',
            expected: {project_name: 'ontologies', model_name: 'uberon', attribute_name: 'name'},
            actual: null
          }],
          missing_ontology_columns: [
            {model: 'labor', template_model: 'labor', column: 'disease'}
          ]
        },
        model_relationships: {
          hierarchy: {
            conforming: false,
            issues: [
              {model: 'tube', template_model: 'sc_seq', parent: 'pool', expected_parent: 'specimen'}
            ]
          },
          single_cell: {
            applicable: true,
            conforming: false,
            child_models: {subject: 'labor', biospecimen: 'specimen', sc_seq: 'tube', sc_seq_pool: 'pool'},
            pool_link: null,
            subject_path: null,
            issues: ['tube has no link attribute to pool'],
            counts: {
              pools: {total: 4, depooled: 3},
              sc_seq: {total: 6, pool_linked: 5, depooled: 3},
              subjects: 2
            }
          }
        }
      }
    });

    global.CONFIG.template_project = 'coprojects_template';

    store = mockStore({
      magma: {models},
      janus: {projects: require('../fixtures/project_names.json')},
      user: {permissions: {labors: {role: 'administrator'}}}
    });

    render(
      <Provider store={store}>
        <ModelMap />
      </Provider>
    );

    fireEvent.click(screen.getByText('Run audit'));

    expect(await screen.findByText('labor → coprojects_template.labor')).toBeTruthy();
    const dialog = within(screen.getByRole('dialog'));
    expect(dialog.getByText('Project audit')).toBeTruthy();
    expect(dialog.queryByText('Needs attention')).toBeNull();
    expect(dialog.queryByText('Conforming')).toBeNull();
    expect(dialog.queryByText('labors')).toBeNull();
    expect(dialog.getByText('Template mappings')).toBeTruthy();
    expect(dialog.getByText('3 mapped')).toBeTruthy();
    expect(dialog.getByText('1 unmapped')).toBeTruthy();
    expect(dialog.queryByText('1 mapped, 1 unmapped, 2 invalid')).toBeNull();
    expect(dialog.getByText('Mapped')).toBeTruthy();
    expect(dialog.getByText('labor → coprojects_template.labor')).toBeTruthy();
    expect(dialog.getByText('Unmapped')).toBeTruthy();
    expect(dialog.getByText('monster')).toBeTruthy();
    expect(dialog.queryByText('labors.monster')).toBeNull();
    expect(dialog.getByText('prize → coprojects_template.award')).toBeTruthy();
    expect(dialog.getByText('cup → template model missing')).toBeTruthy();
    expect(dialog.getByText('Required columns')).toBeTruthy();
    expect(dialog.getByText('labor.name')).toBeTruthy();
    expect(dialog.getByText('labor.species')).toBeTruthy();
    expect(dialog.getByText('Ontology validation')).toBeTruthy();
    expect(dialog.getByText('3 failed')).toBeTruthy();
    expect(dialog.getByText('1 valid')).toBeTruthy();
    expect(dialog.queryByText('human')).toBeNull();
    expect(dialog.getByText('PBMC')).toBeTruthy();
    expect(dialog.getByText('expected ontologies.uberon.name, found none')).toBeTruthy();
    expect(dialog.getByText('column missing')).toBeTruthy();

    fireEvent.click(dialog.getByText('Ontology validation'));
    expect(dialog.queryByText('PBMC')).toBeNull();
    expect(dialog.queryByText('column missing')).toBeNull();
    expect(dialog.getByText('Ontology validation')).toBeTruthy();
    expect(dialog.getByText('3 failed')).toBeTruthy();
    expect(dialog.getByText('1 valid')).toBeTruthy();
    expect(dialog.getByText('labor.name')).toBeTruthy();

    fireEvent.click(dialog.getByText('Ontology validation'));
    expect(dialog.getByText('PBMC')).toBeTruthy();
    expect(dialog.getByText('Single cell models')).toBeTruthy();
    expect(dialog.getByText('Subject relationship: subject → biospecimen → sc_seq')).toBeTruthy();
    expect(dialog.getByText('Pooled relationship: sc_seq → sc_seq_pool')).toBeTruthy();
    expect(dialog.getByText('tube: parent pool, expected specimen')).toBeTruthy();
    expect(dialog.getByText('tube has no link attribute to pool')).toBeTruthy();
    expect(dialog.getByText('2 subjects can be found in search')).toBeTruthy();
    expect(dialog.getByText('2 sc_seq records have a pool but no subject')).toBeTruthy();
    expect(dialog.getByText('1 sc_seq record has no pool')).toBeTruthy();
    expect(dialog.queryByText('Sequencing records')).toBeNull();
    expect(dialog.queryByText('Pools')).toBeNull();
    expect(dialog.queryByText('Subjects reached')).toBeNull();
    expect(dialog.queryByText(/with none/)).toBeNull();
    expect(dialog.queryByText(/depooled|pool_linked/)).toBeNull();
    expect(dialog.queryByText('1/2')).toBeNull();
  });

  it('hides column and ontology checks when the audit has no mapped models', async () => {
    stubUrl({
      verb: 'get',
      path: '/cross_project_search_audit/labors',
      host: 'https://magma.test',
      status: 200,
      response: {
        project: 'labors',
        conforming: false,
        template_mapping: {
          mapped_models: [],
          unmapped_models: ['monster'],
          invalid_mappings: [],
          required_column_results: [
            {model: 'labor', template_model: 'labor', column: 'species', present: false}
          ],
          ontology_value_results: [],
          ontology_valid_count: 1,
          ontology_mismatches: [],
          missing_ontology_columns: []
        },
        model_relationships: {
          hierarchy: {conforming: true, issues: []},
          single_cell: {
            applicable: true,
            conforming: true,
            pool_link: 'pool',
            subject_path: ['tube', 'specimen', 'labor'],
            issues: [],
            counts: {
              pools: {total: 1, depooled: 1},
              sc_seq: {total: 1, pool_linked: 1, depooled: 1},
              subjects: 1
            }
          }
        }
      }
    });

    global.CONFIG.template_project = 'coprojects_template';

    store = mockStore({
      magma: {models},
      janus: {projects: require('../fixtures/project_names.json')},
      user: {permissions: {labors: {role: 'administrator'}}}
    });

    render(
      <Provider store={store}>
        <ModelMap />
      </Provider>
    );

    fireEvent.click(screen.getByText('Run audit'));

    expect(await screen.findByText('No models are mapped. Map models to the template first.')).toBeTruthy();
    const dialog = within(screen.getByRole('dialog'));
    expect(dialog.getByText('Project audit')).toBeTruthy();
    expect(dialog.getByText('1 unmapped')).toBeTruthy();
    expect(dialog.getByText('Mapped')).toBeTruthy();
    expect(dialog.getByText('Unmapped')).toBeTruthy();
    expect(dialog.getByText('monster')).toBeTruthy();
    expect(dialog.queryByText('labors.monster')).toBeNull();
    expect(dialog.queryByText('Required columns')).toBeNull();
    expect(dialog.queryByText('Ontology validation')).toBeNull();
    expect(dialog.queryByText('labor.species')).toBeNull();
    expect(dialog.getByText('Single cell models')).toBeTruthy();
    expect(dialog.getByText('Subject relationship: subject → biospecimen → sc_seq')).toBeTruthy();
    expect(dialog.getByText('Pooled relationship: sc_seq → sc_seq_pool')).toBeTruthy();
    expect(dialog.getByText('1 subject can be found in search')).toBeTruthy();
    expect(dialog.queryByText(/sc_seq record/)).toBeNull();
    expect(dialog.queryByText('Sequencing records')).toBeNull();
    expect(dialog.queryByText('Pools')).toBeNull();
    expect(dialog.queryByText('Subjects reached')).toBeNull();
    expect(dialog.queryByText(/with none/)).toBeNull();
    expect(dialog.queryByText(/linked to a pool only/)).toBeNull();
    expect(dialog.queryByText(/not linked to a pool/)).toBeNull();
    expect(dialog.queryByText(/depooled|pool_linked/)).toBeNull();
  });

  it('hides the project audit button on the template project', async () => {
    global.CONFIG.project_name = 'coprojects_template';
    global.CONFIG.template_project = 'coprojects_template';

    store = mockStore({
      magma: {models},
      janus: {projects: require('../fixtures/project_names.json')},
      user: {permissions: {coprojects_template: {role: 'administrator'}}}
    });

    render(
      <Provider store={store}>
        <ModelMap />
      </Provider>
    );

    await waitFor(() => screen.findByText('Attribute'));
    expect(screen.queryByText('Run audit')).toBeNull();
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
