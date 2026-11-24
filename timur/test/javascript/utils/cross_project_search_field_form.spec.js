import {
  TYPE_HELP,
  emptyForm,
  emptyModalityForm,
  fieldGroups,
  formBody,
  formError,
  formFromRow,
  modalityFormBody,
  modalityFormFromRow,
  newFieldForm,
  normalizeFields,
  normalizeModalities,
  countPill,
  parseQueryText,
  presenceQuery,
  queryText,
  queryTextError,
  showUnderOptions
} from '../../../lib/client/jsx/components/cross_project_search/field_form';

const row = {id: 7, label: 'Species', query: ['subject', '::distinct', 'species'], kind: 'filter', modality_id: null, modality: null, position: 2};
const count = {id: 8, label: 'Subjects', query: ['subject', '::count'], kind: 'count', modality_id: null, modality: null, position: 1};
const pools = {id: 12, label: 'Pools', query: ['pool_model', '::count'], kind: 'count', modality_id: 9, modality: 'single cell', position: 6};
const cell = {id: 9, name: 'single cell', model: 'assay_model', position: 1};
const bulk = {id: 10, name: 'bulk RNA-seq', model: 'other_model', position: 2};
const modalities = [bulk, cell];

describe('normalizeFields / normalizeModalities', () => {
  it('turns field rows into search fields', () => {
    expect(normalizeFields([row, count, pools]).map((f) => f.name)).toEqual(['7', '8', '12']);
    expect(normalizeFields([row])[0]).toEqual({name: '7', label: 'Species', query: row.query, kind: 'filter', modalityId: null});
    expect(normalizeFields([pools])[0]).toMatchObject({kind: 'count', modalityId: 9});
    expect(normalizeFields(null)).toEqual([]);
  });

  it('derives the presence query from the model string', () => {
    expect(presenceQuery('anything')).toEqual(['anything', '::count']);
    expect(normalizeModalities(modalities)).toEqual([
      {id: 9, name: 'single cell', model: 'assay_model', position: 1, query: ['assay_model', '::count']},
      {id: 10, name: 'bulk RNA-seq', model: 'other_model', position: 2, query: ['other_model', '::count']}
    ]);
    expect(normalizeModalities(undefined)).toEqual([]);
  });
});

describe('query text', () => {
  it('is compact JSON and parses back; invalid text is passed through raw', () => {
    expect(queryText(row.query)).toEqual('["subject","::distinct","species"]');
    expect(parseQueryText('["subject", "::count"]')).toEqual(['subject', '::count']);
    expect(parseQueryText('subject')).toEqual('subject');
    expect(parseQueryText('{"a":1}')).toEqual({a: 1});
  });

  it('reports invalid JSON for the form', () => {
    expect(queryTextError('["a"]')).toBeNull();
    expect(queryTextError('[')).toMatch(/not valid JSON/);
    expect(formError({queryText: '[1'})).toMatch(/not valid JSON/);
    expect(formError({queryText: '["a"]'})).toBeNull();
  });
});

describe('form', () => {
  it('round trips a field and always sends the kind', () => {
    expect(formBody(formFromRow(row))).toEqual({kind: 'filter', label: 'Species', query: row.query, modality_id: null, position: 2});
    expect(formBody(formFromRow(pools))).toEqual({kind: 'count', label: 'Pools', query: pools.query, modality_id: 9, position: 6});
    expect(formFromRow(pools).modalityId).toEqual('9');
  });

  it('a new form is a count, Global, with no position; Show under can be preset', () => {
    expect(emptyForm()).toMatchObject({id: null, kind: 'count', modalityId: '', position: ''});
    const form = {...newFieldForm(9), label: 'x', queryText: '["a","::count"]'};
    expect(formBody(form)).toEqual({kind: 'count', label: 'x', query: ['a', '::count'], modality_id: 9});
    expect(newFieldForm(null).modalityId).toEqual('');
  });

  it('sends raw text for invalid JSON and a non-numeric position as text', () => {
    const body = formBody({...emptyForm(), label: 'x', queryText: 'oops', position: 'two'});
    expect(body).toMatchObject({query: 'oops', position: 'two', modality_id: null});
  });

  it('has one line of helper text per type, naming no model', () => {
    expect(TYPE_HELP.count).toEqual('Must end with ::count. Example: ["<model>","::count"]');
    expect(TYPE_HELP.filter).toEqual('Must end with ::distinct or ::first. Example: ["<model>","::distinct","<attribute>"]');
  });

  it('round trips a modality form as name and model', () => {
    expect(emptyModalityForm()).toEqual({id: null, name: '', model: ''});
    expect(modalityFormBody(modalityFormFromRow(cell))).toEqual({name: 'single cell', model: 'assay_model'});
  });
});

describe('showUnderOptions', () => {
  it('is Global, then each modality name in position order', () => {
    expect(showUnderOptions(modalities)).toEqual([
      {value: '', label: 'Global'},
      {value: '9', label: 'single cell'},
      {value: '10', label: 'bulk RNA-seq'}
    ]);
    expect(showUnderOptions(null)).toEqual([{value: '', label: 'Global'}]);
  });
});

describe('fieldGroups', () => {
  const plates = {id: 13, label: 'Plates', query: ['plate_model', '::count'], kind: 'count', modality_id: 10, position: 4};
  const stray = {id: 14, label: 'Orphan', query: ['x', '::count'], kind: 'count', modality_id: 99, position: 3};

  it('is Global, then one group per modality, fields in position order', () => {
    const groups = fieldGroups([row, count, pools, plates, stray], modalities);
    expect(groups.map((group) => ({title: group.title, count: group.count, labels: group.rows.map((item) => item.label)}))).toEqual([
      {title: 'Global', count: 3, labels: ['Subjects', 'Species', 'Orphan']},
      {title: 'single cell', count: 1, labels: ['Pools']},
      {title: 'bulk RNA-seq', count: 1, labels: ['Plates']}
    ]);
  });

  it('still has a Global group when there is nothing', () => {
    expect(fieldGroups([], []).map((group) => group.title)).toEqual(['Global']);
    expect(fieldGroups(null, null)[0].count).toEqual(0);
  });
});

describe('countPill', () => {
  it('pluralizes', () => {
    expect(countPill(1, 'field')).toEqual('1 field');
    expect(countPill(9, 'field')).toEqual('9 fields');
  });
});
