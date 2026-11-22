import React from 'react';
import {Provider} from 'react-redux';
import {render, screen, fireEvent, within} from '@testing-library/react';
import '@testing-library/jest-dom/extend-expect';
import {cleanStubs, mockStore, stubUrl} from '../helpers';
import CrossProjectSearchModal from '../../../lib/client/jsx/components/cross_project_search/search_modal';

const host = 'https://magma.test';

const fieldRows = [
  {id: 1, label: 'Species', query: ['model_a', '::distinct', 'species'], kind: 'filter', modality_id: null, modality: null, position: 1},
  {id: 2, label: 'Records', query: ['model_a', '::count'], kind: 'count', modality_id: null, modality: null, position: 2}
];
const modalityRows = [{id: 5, name: 'Assay one', model: 'model_b', position: 1}];
const queries = [['model_a', '::distinct', 'species'], ['model_a', '::count'], ['model_b', '::count']];

const stubRows = () => {
  stubUrl({verb: 'get', path: '/cross_project_search/fields', host, response: {fields: fieldRows}});
  stubUrl({verb: 'get', path: '/cross_project_search/modalities', host, response: {modalities: modalityRows}});
};

const renderModal = () => render(
  <Provider store={mockStore({})}>
    <CrossProjectSearchModal open={true} onClose={() => {}} />
  </Provider>
);

describe('CrossProjectSearchModal', () => {
  beforeEach(() => {
    global.CONFIG = {magma_host: host};
    window.localStorage.clear();
  });

  afterEach(cleanStubs);

  it('loads fields and modalities, runs one query batch and shows a row per project', async () => {
    stubRows();
    stubUrl({
      verb: 'post',
      path: '/cross_project_search/query',
      host,
      request: {queries},
      response: {projects: {zeta: [['Mus musculus', 'Homo sapiens'], 3, 2], alpha: [null, null, null]}}
    });

    renderModal();

    expect(await screen.findByText('zeta')).toBeTruthy();
    expect(screen.getByText('2 of 2 projects')).toBeTruthy();
    const zeta = within(screen.getByText('zeta').closest('tr'));
    expect(zeta.getByText('Homo sapiens, Mus musculus')).toBeTruthy();
    expect(zeta.getByText('Assay one')).toBeTruthy();
    expect(zeta.getByText('3')).toBeTruthy();
    expect(within(screen.getByText('alpha').closest('tr')).getAllByText('n/a')).toHaveLength(3);
  });

  it('shows the server error from a failed batch and keeps the Fields view usable', async () => {
    stubRows();
    stubUrl({
      verb: 'post',
      path: '/cross_project_search/query',
      host,
      request: {queries},
      status: 422,
      response: {error: "Query 3: Unknown model 'model_b'"}
    });

    renderModal();

    expect(await screen.findByText("Query 3: Unknown model 'model_b'")).toBeTruthy();
    fireEvent.click(screen.getByRole('tab', {name: 'Fields'}));
    expect(screen.getByRole('tab', {name: 'Modalities'})).toHaveAttribute('aria-selected', 'true');
    expect(screen.getByText('Assay one')).toBeTruthy();
    expect(screen.getByText('model_b')).toBeTruthy();

    const [, fieldsTab] = screen.getAllByRole('tab', {name: 'Fields'});
    fireEvent.click(fieldsTab);
    expect(screen.getByText('Species')).toBeTruthy();
    expect(screen.getByText('["model_a","::count"]')).toBeTruthy();
  });
});
