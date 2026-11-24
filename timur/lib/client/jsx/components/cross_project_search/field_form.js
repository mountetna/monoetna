// Pure helpers for the Fields admin view and for turning server rows into
// search fields. No React, no fetching.
//
// Modalities come from GET /cross_project_search/modalities:
//   {id, name, model, position, updated_at}
// A project has a modality when a count of `model` is a number > 0. That
// presence query is derived here and is never stored.
//
// Fields come from GET /cross_project_search/fields:
//   {id, label, query, kind, modality_id, modality, position, updated_at}
// kind is 'count' or 'filter'. modality_id null is global.

import {COUNT, FILTER} from './facet_filter';

const byPosition = (items) => (
  items.map((item, index) => ({item, index}))
    .sort((a, b) => {
      const pa = typeof a.item.position === 'number' ? a.item.position : null;
      const pb = typeof b.item.position === 'number' ? b.item.position : null;
      if (pa !== null && pb !== null && pa !== pb) return pa - pb;
      return a.index - b.index;
    })
    .map((entry) => entry.item)
);

// ["<model>", "::count"] from whatever model string the row has.
export const presenceQuery = (model) => [model, '::count'];

// Server field rows become the {name, label, query, kind, modalityId} fields
// the search logic uses.
export const normalizeFields = (rows) => (
  (rows || []).map((row) => ({
    name: String(row.id),
    label: row.label,
    query: row.query,
    kind: row.kind === COUNT ? COUNT : FILTER,
    modalityId: row.modality_id == null ? null : row.modality_id
  }))
);

// Server modality rows become {id, name, model, position, query}. `query` is
// the derived presence query, sent in the same batch as the field queries.
export const normalizeModalities = (rows) => (
  byPosition(rows || []).map((row) => ({
    id: row.id,
    name: row.name,
    model: row.model,
    position: row.position,
    query: presenceQuery(row.model)
  }))
);

// The Fields list, in position order: Global, then one group per modality.
// A field whose modality is missing is shown under Global so it can be edited.
export const fieldGroups = (fieldRows, modalities) => {
  const fields = fieldRows || [];
  const mods = byPosition(modalities || []);
  const ids = mods.map((modality) => modality.id);
  const modalityOf = (row) => (
    row.modality_id != null && ids.includes(row.modality_id) ? row.modality_id : null
  );
  const under = (id) => byPosition(fields.filter((row) => modalityOf(row) === id));
  return [{key: 'global', title: 'Global', rows: under(null)}]
    .concat(mods.map((modality) => ({
      key: `modality-${modality.id}`,
      title: modality.name,
      rows: under(modality.id)
    })))
    .map((group) => ({...group, count: group.rows.length}));
};

// Compact JSON text of a query for lists and the edit form.
export const queryText = (query) => JSON.stringify(query);

// The text area holds JSON. Valid JSON is sent as the parsed value (the server
// checks it is a proper query); anything else is sent as the raw text so the
// server answers with its own message.
export const parseQueryText = (text) => {
  try {
    return JSON.parse(text);
  } catch (e) {
    return text;
  }
};

// null when the text is valid JSON, otherwise a message for the form.
export const queryTextError = (text) => {
  try {
    JSON.parse(text);
    return null;
  } catch (e) {
    return `Query is not valid JSON: ${e.message}`;
  }
};

export const TYPES = [
  {value: COUNT, label: 'Count'},
  {value: FILTER, label: 'Filter'}
];

export const emptyForm = () => ({id: null, kind: COUNT, label: '', queryText: '[]', modalityId: '', position: ''});

// A new field form. modalityId null/'' is Global.
export const newFieldForm = (modalityId) => ({
  ...emptyForm(),
  modalityId: modalityId == null ? '' : String(modalityId)
});

export const TYPE_HELP = {
  [COUNT]: 'Must end with ::count. Example: ["<model>","::count"]',
  [FILTER]: 'Must end with ::distinct or ::first. Example: ["<model>","::distinct","<attribute>"]'
};

export const formFromRow = (row) => ({
  id: row.id,
  kind: row.kind,
  label: row.label,
  queryText: queryText(row.query),
  modalityId: row.modality_id == null ? '' : String(row.modality_id),
  position: String(row.position)
});

// The body of POST /cross_project_search/fields and /cross_project_search/fields/:id. The kind is always sent. An empty
// position is left out (new rows go to the end; edits keep their place).
export const formBody = (form) => {
  const body = {
    kind: form.kind,
    label: form.label,
    query: parseQueryText(form.queryText),
    modality_id: String(form.modalityId || '') === '' ? null : Number(form.modalityId)
  };
  const position = String(form.position).trim();
  if (position !== '') body.position = /^-?\d+$/.test(position) ? Number(position) : position;
  return body;
};

export const emptyModalityForm = () => ({id: null, name: '', model: ''});

export const modalityFormFromRow = (row) => ({
  id: row.id,
  name: row.name,
  model: row.model
});

// The body of POST /cross_project_search/modalities and /cross_project_search/modalities/:id.
export const modalityFormBody = (form) => ({
  name: form.name,
  model: form.model
});

// Options of the Show under select: Global, then every modality by name.
export const showUnderOptions = (modalities) => (
  [{value: '', label: 'Global'}].concat(
    byPosition(modalities || []).map((modality) => ({value: String(modality.id), label: modality.name}))
  )
);

// "1 field", "9 fields" — the small pill beside a group header.
export const countPill = (count, noun, plural) => `${count} ${count === 1 ? noun : (plural || `${noun}s`)}`;

// Local check before sending: valid JSON. null when fine.
export const formError = (form) => queryTextError(form.queryText);

