// A query batch is every field query in field order, then every modality's
// presence query in modality order; each project's results follow that order.
// A field's kind, never the shape of its answer, decides how it is drawn.

const isCount = (value) => typeof value === 'number';

export const COUNT = 'count';
export const FILTER = 'filter';
export const PROJECT = 'project';

const modalityName = (modalities, id) => {
  const modality = (modalities || []).find((item) => item.id === id);
  return modality ? modality.name : null;
};

// Global fields are always in scope. A field with a modality is in scope only
// when that modality's name is among `activeNames`.
const inScope = (field, modalities, activeNames) => {
  if (field.modalityId == null) return true;
  const name = modalityName(modalities, field.modalityId);
  return name != null && (activeNames || []).includes(name);
};

export const countFields = (fields) => fields.filter((field) => field.kind === COUNT);

// The derived Modality facet/column.
export const MODALITY_ID = 'modality';
export const MODALITY_LABEL = 'Modality';

// Names of the modalities whose answer is a number > 0, in modality order.
// answers[i] belongs to modalities[i]. null is not a match.
export const projectModalities = (modalities, answers) => (
  (modalities || [])
    .filter((modality, index) => isCount(answers[index]) && answers[index] > 0)
    .map((modality) => modality.name)
);

// Filter fields nested under one modality, in the given order.
export const modalityFilterFields = (fields, modality) => (
  (fields || []).filter((field) => field.kind !== COUNT && field.modalityId != null && field.modalityId === modality.id)
);

// The Modality facet's selection keys: the modality names plus every nested
// filter field. Everything ticked under these keys is one OR group.
export const modalityGroup = (fields, modalities) => (
  (modalities || []).length
    ? [MODALITY_ID].concat(...modalities.map((modality) => modalityFilterFields(fields, modality).map((field) => field.name)))
    : []
);

// Names of the modalities in play: ticked themselves or with a nested value ticked.
export const activeModalities = (selected, fields, modalities) => {
  const current = selected || {};
  const ticked = (name) => (current[name] || []).length > 0;
  return (modalities || [])
    .filter((modality) => (current[MODALITY_ID] || []).includes(modality.name)
      || modalityFilterFields(fields, modality).some((field) => ticked(field.name)))
    .map((modality) => modality.name);
};

// The Modality facet as a field-like object, first in the sidebar; [] when there are no modalities.
export const modalityFacet = (modalities) => (
  (modalities || []).length ? [{name: MODALITY_ID, label: MODALITY_LABEL, kind: FILTER}] : []
);

// Counts shown now: global ones always; modality ones when their modality is active.
export const applicableCounts = (fields, modalities, activeNames = []) => (
  countFields(fields).filter((field) => inScope(field, modalities, activeNames))
);

export const facetFields = (fields) => fields.filter((field) => field.kind !== COUNT);

// Filter columns shown now: global ones always; modality ones when their
// modality is active (see activeModalities).
export const applicableFilters = (fields, modalities, activeNames = []) => (
  facetFields(fields).filter((field) => inScope(field, modalities, activeNames))
);

// The sidebar: Modality first (with its nested filters), then the global filters.
export const sidebarFacets = (fields, modalities) => (
  modalityFacet(modalities).concat(facetFields(fields).filter((field) => field.modalityId == null))
);

const displayValue = (value) => (typeof value === 'string' ? value.trim() : String(value));

const compareValues = (left, right) => (
  left.localeCompare(right, undefined, {numeric: true, sensitivity: 'base'})
);

// The distinct values of a project for one field as sorted strings. Booleans
// become 'true'/'false'; null and empty values are dropped; a number is never a value.
export const projectValues = (project, name) => {
  const raw = project.values ? project.values[name] : null;
  const list = Array.isArray(raw) ? raw : [raw];
  const texts = list
    .filter((item) => item != null && !isCount(item))
    .map(displayValue)
    .filter((text) => text !== '');
  return Array.from(new Set(texts)).sort(compareValues);
};

// Distinct values of one facet across all projects, each with the number of
// projects that have it.
export const facetOptions = (name, projects) => {
  const counts = {};
  projects.forEach((project) => {
    projectValues(project, name).forEach((value) => {
      counts[value] = (counts[value] || 0) + 1;
    });
  });
  return Object.keys(counts)
    .sort(compareValues)
    .map((value) => ({value, count: counts[value]}));
};

// What the Modality column shows, most specific first: for each modality in
// order, the project's values for that modality's nested filters when it has
// any, else the modality's name when the project has the modality. Sorted
// within a modality, duplicates dropped.
export const modalityDetail = (project, fields, modalities) => {
  const has = projectValues(project, MODALITY_ID);
  const parts = [];
  (modalities || []).forEach((modality) => {
    const nested = Array.from(new Set(
      [].concat(...modalityFilterFields(fields, modality).map((field) => projectValues(project, field.name)))
    )).sort(compareValues);
    if (nested.length) parts.push(...nested);
    else if (has.includes(modality.name)) parts.push(modality.name);
  });
  return Array.from(new Set(parts));
};

// One row per project: {project_name, values: {fieldName: result, modality: [names]}, modalityDetail}.
// Each project's result list is the fields' answers followed by the modalities' answers.
export const buildProjects = (fields, modalities, results) => (
  Object.keys(results || {}).sort().map((projectName) => {
    const row = results[projectName] || [];
    const values = {};
    fields.forEach((field, index) => {
      values[field.name] = row[index] === undefined ? null : row[index];
    });
    const project = {project_name: projectName, values};
    if ((modalities || []).length) {
      values[MODALITY_ID] = projectModalities(
        modalities,
        modalities.map((modality, index) => row[fields.length + index])
      );
      project.modalityDetail = modalityDetail(project, fields, modalities);
    }
    return project;
  })
);

const matchesFacet = (project, name, values) => {
  const own = projectValues(project, name);
  return values.some((value) => own.includes(value));
};

// selected: {name: [values]}. AND across facets, OR within a facet. Names in
// `group` together act as one facet (OR across all their ticked values). A
// project matches a ticked value if any of its values equals it. A facet with
// nothing ticked does not filter.
export const filterProjects = (projects, selected, group = []) => {
  const current = selected || {};
  const active = Object.keys(current).filter((name) => current[name].length > 0);
  const andNames = active.filter((name) => !group.includes(name));
  const orNames = active.filter((name) => group.includes(name));
  return projects.filter((project) => (
    andNames.every((name) => matchesFacet(project, name, current[name]))
    && (orNames.length === 0 || orNames.some((name) => matchesFacet(project, name, current[name])))
  ));
};

export const summaryText = (shown, total) => `${shown.length} of ${total.length} projects`;

export const toggleValue = (selected, name, value) => {
  const current = selected[name] || [];
  const next = current.includes(value)
    ? current.filter((item) => item !== value)
    : current.concat([value]);
  const copy = {...selected};
  if (next.length) copy[name] = next;
  else delete copy[name];
  return copy;
};

// Project, Modality (when there are modalities), the filters that apply, then the counts
// that apply to the active modalities (see activeModalities), each in the position order given.
export const tableColumns = (fields, modalities, activeNames = []) => (
  [{id: PROJECT, label: 'Project', kind: PROJECT}]
    .concat(modalityFacet(modalities).map((field) => ({id: field.name, label: field.label, kind: FILTER})))
    .concat(applicableFilters(fields, modalities, activeNames).map((field) => ({id: field.name, label: field.label, kind: FILTER})))
    .concat(applicableCounts(fields, modalities, activeNames).map((field) => ({id: field.name, label: field.label, kind: COUNT})))
);

// Text of one cell. Counts are numbers, filters are their distinct values
// joined by `separator`; anything missing or of the wrong shape is 'n/a'.
export const cellValue = (project, column, separator = ', ') => {
  if (column.kind === PROJECT) return project.project_name;
  if (column.kind === FILTER) {
    const values = column.id === MODALITY_ID && project.modalityDetail
      ? project.modalityDetail
      : projectValues(project, column.id);
    return values.length ? values.join(separator) : 'n/a';
  }
  const value = project.values[column.id];
  return isCount(value) ? String(value) : 'n/a';
};

export const tableRows = (projects, columns, separator) => (
  projects.map((project) => columns.map((column) => cellValue(project, column, separator)))
);

// Same cells as tableRows, but list values are joined with "; " so the
// commas inside a value are not confused with the CSV's own.
export const csvRows = (projects, columns) => tableRows(projects, columns, '; ');

const csvEscape = (value) => {
  const text = value == null ? '' : String(value);
  if (/[",\n\r]/.test(text)) return `"${text.replace(/"/g, '""')}"`;
  return text;
};

// Exactly the rows and columns on screen. Header is the column labels.
export const buildCsv = (columns, rows) => (
  [columns.map((column) => column.label)]
    .concat(rows)
    .map((row) => row.map(csvEscape).join(','))
    .join('\n')
);

// Live facet options: every value that exists for the facet, each with the
// number of projects that would match if it were ticked, given the OTHER
// facets' selections (the facet's own selection is ignored, so siblings stay
// available for OR). Count 0 values are kept so the UI can grey them out.
// `group` is the Modality facet's OR group (see modalityGroup).
export const liveFacetOptions = (name, projects, selected, group = []) => {
  const others = {...(selected || {})};
  delete others[name];
  const matching = filterProjects(projects, others, group);
  const live = {};
  facetOptions(name, matching).forEach(({value, count}) => { live[value] = count; });
  return facetOptions(name, projects).map(({value}) => ({value, count: live[value] || 0}));
};

// The Modality facet: every modality in the given order with its live count,
// each with `groups`: one per nested filter field, {name, label, options}, in
// field order, options being every value of that field with its live count.
// Counts ignore everything ticked in the Modality facet itself.
export const liveModalityOptions = (modalities, projects, selected, fields = []) => {
  const group = modalityGroup(fields, modalities);
  const others = {...(selected || {})};
  group.forEach((name) => { delete others[name]; });
  const matching = filterProjects(projects, others);
  const liveCounts = (name) => {
    const live = {};
    facetOptions(name, matching).forEach(({value, count}) => { live[value] = count; });
    return live;
  };
  const liveModality = liveCounts(MODALITY_ID);
  return (modalities || []).map((modality) => ({
    value: modality.name,
    count: liveModality[modality.name] || 0,
    groups: modalityFilterFields(fields, modality).map((field) => {
      const live = liveCounts(field.name);
      return {
        name: field.name,
        label: field.label,
        options: facetOptions(field.name, projects).map(({value}) => ({value, count: live[value] || 0}))
      };
    })
  }));
};

export const formatNumber = (number) => number.toLocaleString('en-US');

// Sum of a count column over the given projects; null/n/a ignored.
// Returns null when no project has a number.
export const columnTotal = (projects, column) => {
  const numbers = projects.map((project) => project.values[column.id]).filter(isCount);
  return numbers.length ? numbers.reduce((sum, number) => sum + number, 0) : null;
};

// Cells for the bottom row: 'Total', then each count column's sum, blank for others.
export const totalsRow = (projects, columns) => (
  columns.map((column, index) => {
    if (index === 0) return 'Total';
    if (column.kind !== COUNT) return '';
    const total = columnTotal(projects, column);
    return total === null ? '' : formatNumber(total);
  })
);

// ---- Column visibility. Order is never user-set: it always follows Position.

export const HIDDEN_COLUMNS_KEY = 'timur/crossProjectSearch/hiddenColumns';

// The stored hidden-column text, or null.
export const readHiddenColumns = (storage) => {
  try {
    return (storage || window.localStorage).getItem(HIDDEN_COLUMNS_KEY);
  } catch (e) {
    return null;
  }
};

export const writeHiddenColumns = (text, storage) => {
  try {
    (storage || window.localStorage).setItem(HIDDEN_COLUMNS_KEY, text);
  } catch (e) {
    // Storage may be unavailable (private mode, quota); the choice then lasts for the session only.
  }
};

// Columns the user may toggle: every field column, never Project.
export const toggleableColumns = (columns) => columns.filter((column) => column.kind !== PROJECT);

// Stored text -> list of hidden ids that still exist. Anything unreadable is [].
export const parseHiddenColumns = (text, fields, modalities) => {
  let parsed;
  try { parsed = JSON.parse(text); } catch (error) { return []; }
  if (!Array.isArray(parsed)) return [];
  const known = (fields || []).map((field) => field.name).concat((modalities || []).length ? [MODALITY_ID] : []);
  return Array.from(new Set(parsed.map(String))).filter((id) => known.includes(id));
};

export const toggleHidden = (hidden, id) => (
  hidden.includes(id) ? hidden.filter((item) => item !== id) : hidden.concat([id])
);

// The columns to draw, export and total: Project always, then all not hidden.
export const visibleColumns = (columns, hidden) => (
  columns.filter((column) => column.kind === PROJECT || !(hidden || []).includes(column.id))
);
