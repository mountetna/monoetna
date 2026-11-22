import {
  buildCsv,
  csvRows,
  liveFacetOptions,
  totalsRow,
  buildProjects,
  countFields,
  facetFields,
  facetOptions,
  filterProjects,
  projectValues,
  summaryText,
  tableColumns,
  projectModalities,
  modalityFacet,
  liveModalityOptions,
  applicableFilters,
  activeModalities,
  modalityDetail,
  modalityFilterFields,
  modalityGroup,
  sidebarFacets,
  applicableCounts,
  tableRows,
  toggleValue,
  parseHiddenColumns,
  readHiddenColumns,
  writeHiddenColumns,
  HIDDEN_COLUMNS_KEY,
  toggleHidden,
  toggleableColumns,
  visibleColumns
} from '../../../lib/client/jsx/components/cross_project_search/facet_filter';

// Search fields as normalizeFields builds them.
const fields = [
  {name: 'subjects', label: 'Subjects', query: ['subject', '::count'], kind: 'count'},
  {name: 'species', label: 'Species', query: ['subject', '::distinct', 'species'], kind: 'filter'},
  {name: 'tissue', label: 'Tissue', query: ['sample', '::distinct', 'tissue'], kind: 'filter'},
  {name: 'completed', label: 'Completed', query: ['project', '::first', 'completed'], kind: 'filter'},
  {name: 'funding', label: 'Funding Source', query: ['project', '::first', 'funding_source'], kind: 'filter'},
  {name: 'samples', label: 'Samples', query: ['sample', '::count'], kind: 'count'}
];

// As POST /cross_project_search/query sends them: one result per field, in order.
const results = {
  example_project: [3, ['Homo sapiens', 'Mus musculus'], ['blood'], false, 'NIH, "R01"', 10],
  zebra_atlas: [7, ['Danio rerio'], ['blood', 'gill'], true, null, 4],
  ontologies: [null, null, null, null, null, null]
};

const projects = buildProjects(fields, [], results);
const names = (rows) => rows.map((project) => project.project_name);
const labels = (columns) => columns.map((column) => column.label);

describe('buildProjects', () => {
  it('makes one sorted row per project keyed by field name', () => {
    expect(names(projects)).toEqual(['example_project', 'ontologies', 'zebra_atlas']);
    expect(projects[0].values.species).toEqual(['Homo sapiens', 'Mus musculus']);
    expect(projects[0].values.completed).toBe(false);
  });

  it('treats a missing result as null', () => {
    expect(buildProjects(fields, [], {short: [1]})[0].values.species).toBeNull();
  });
});

describe('facetFields and countFields', () => {
  it('split fields by kind', () => {
    expect(facetFields(fields).map((field) => field.name)).toEqual(['species', 'tissue', 'completed', 'funding']);
    expect(countFields(fields).map((field) => field.name)).toEqual(['subjects', 'samples']);
  });
});

describe('projectValues', () => {
  it('lists strings, booleans as true/false and list members, sorted; never null or numbers', () => {
    expect(projectValues(projects[0], 'species')).toEqual(['Homo sapiens', 'Mus musculus']);
    expect(projectValues({values: {x: ['b', 'A', 'c']}}, 'x')).toEqual(['A', 'b', 'c']);
    expect(projectValues(projects[0], 'completed')).toEqual(['false']);
    expect(projectValues(projects[2], 'completed')).toEqual(['true']);
    expect(projectValues(projects[2], 'funding')).toEqual([]);
    expect(projectValues(projects[0], 'subjects')).toEqual([]);
    expect(projectValues({values: {x: [null, '', 'a', 'a']}}, 'x')).toEqual(['a']);
  });
});

describe('facetOptions', () => {
  it('counts projects per distinct value and never lists null', () => {
    expect(facetOptions('tissue', projects)).toEqual([
      {value: 'blood', count: 2},
      {value: 'gill', count: 1}
    ]);
    expect(facetOptions('completed', projects)).toEqual([
      {value: 'false', count: 1},
      {value: 'true', count: 1}
    ]);
    expect(facetOptions('funding', projects)).toEqual([{value: 'NIH, "R01"', count: 1}]);
  });

  it('gives a number field no options', () => {
    expect(facetOptions('subjects', projects)).toEqual([]);
  });
});

describe('filterProjects', () => {
  it('returns everything when nothing is ticked', () => {
    expect(filterProjects(projects, {})).toHaveLength(3);
    expect(filterProjects(projects, {species: []})).toHaveLength(3);
  });

  it('matches a project when any of its values equals the ticked value', () => {
    expect(names(filterProjects(projects, {species: ['Homo sapiens']}))).toEqual(['example_project']);
    expect(names(filterProjects(projects, {species: ['Mus musculus']}))).toEqual(['example_project']);
  });

  it('is OR within a facet', () => {
    expect(names(filterProjects(projects, {species: ['Homo sapiens', 'Danio rerio']})))
      .toEqual(['example_project', 'zebra_atlas']);
  });

  it('is AND across facets', () => {
    const selected = {species: ['Homo sapiens', 'Danio rerio'], tissue: ['gill']};
    expect(names(filterProjects(projects, selected))).toEqual(['zebra_atlas']);
    expect(filterProjects(projects, {species: ['Homo sapiens'], tissue: ['gill']})).toEqual([]);
  });

  it('filters on a boolean shown as true/false', () => {
    expect(names(filterProjects(projects, {completed: ['true']}))).toEqual(['zebra_atlas']);
    expect(names(filterProjects(projects, {completed: ['false']}))).toEqual(['example_project']);
  });

  it('drops projects with no value for a ticked facet', () => {
    expect(names(filterProjects(projects, {funding: ['NIH, "R01"']}))).toEqual(['example_project']);
  });
});

describe('summaryText and toggleValue', () => {
  it('reads N of M projects', () => {
    expect(summaryText(filterProjects(projects, {species: ['Danio rerio']}), projects)).toEqual('1 of 3 projects');
  });

  it('adds and removes values and drops empty facets', () => {
    let selected = toggleValue({}, 'species', 'Homo sapiens');
    selected = toggleValue(selected, 'species', 'Danio rerio');
    expect(selected).toEqual({species: ['Homo sapiens', 'Danio rerio']});
    selected = toggleValue(selected, 'species', 'Homo sapiens');
    selected = toggleValue(selected, 'species', 'Danio rerio');
    expect(selected).toEqual({});
  });
});

describe('columns', () => {
  it('is Project, every filter field, then every count field, in order', () => {
    const columns = tableColumns(fields, []);
    expect(labels(columns)).toEqual(['Project', 'Species', 'Tissue', 'Completed', 'Funding Source', 'Subjects', 'Samples']);
    expect(tableRows(projects, columns)).toEqual([
      ['example_project', 'Homo sapiens, Mus musculus', 'blood', 'false', 'NIH, "R01"', '3', '10'],
      ['ontologies', 'n/a', 'n/a', 'n/a', 'n/a', 'n/a', 'n/a'],
      ['zebra_atlas', 'Danio rerio', 'blood, gill', 'true', 'n/a', '7', '4']
    ]);
  });

  it('has only Project when there are no fields', () => {
    expect(labels(tableColumns([], []))).toEqual(['Project']);
  });
});

describe('kind mismatch', () => {
  const mismatched = [
    {name: 'c', label: 'C', query: ['x', '::count'], kind: 'count'},
    {name: 'f', label: 'F', query: ['x', '::distinct', 'a'], kind: 'filter'}
  ];
  const odd = buildProjects(mismatched, [], {p: [['a', 'b'], 5], q: ['text', ['v']]});

  it('shows n/a for a count that is not a number or a filter that is a number', () => {
    expect(tableRows(odd, tableColumns(mismatched, []))).toEqual([['p', 'n/a', 'n/a'], ['q', 'v', 'n/a']]);
    expect(facetOptions('f', odd)).toEqual([{value: 'v', count: 1}]);
  });
});

describe('csv', () => {
  it('exports exactly the visible rows and columns, joining lists with "; " and quoting', () => {
    const columns = tableColumns(fields, []);
    const shown = filterProjects(projects, {tissue: ['blood']});
    expect(buildCsv(columns, csvRows(shown, columns))).toEqual([
      'Project,Species,Tissue,Completed,Funding Source,Subjects,Samples',
      'example_project,Homo sapiens; Mus musculus,blood,false,"NIH, ""R01""",3,10',
      'zebra_atlas,Danio rerio,blood; gill,true,n/a,7,4'
    ].join('\n'));
    expect(csvRows(projects, columns)).toHaveLength(3);
    expect(buildCsv([{label: 'A "b"'}], [['x\ny']])).toEqual('"A ""b"""\n"x\ny"');
  });
});

describe('liveFacetOptions', () => {
  it('counts against the other facets and ignores its own selection', () => {
    expect(liveFacetOptions('tissue', projects, {})).toEqual([
      {value: 'blood', count: 2},
      {value: 'gill', count: 1}
    ]);
    expect(liveFacetOptions('tissue', projects, {species: ['Homo sapiens']})).toEqual([
      {value: 'blood', count: 1},
      {value: 'gill', count: 0}
    ]);
    // own selection ignored: siblings stay counted for OR
    expect(liveFacetOptions('tissue', projects, {tissue: ['gill']})).toEqual([
      {value: 'blood', count: 2},
      {value: 'gill', count: 1}
    ]);
  });

  it('is AND across the other facets', () => {
    const selected = {species: ['Danio rerio'], completed: ['false']};
    expect(liveFacetOptions('tissue', projects, selected)).toEqual([
      {value: 'blood', count: 0},
      {value: 'gill', count: 0}
    ]);
  });
});

describe('totals row', () => {
  const columns = tableColumns(fields, []);

  it('formats thousands', () => {
    const big = buildProjects(fields, [], {a: [2000, null, null, null, null, null], b: [1420, null, null, null, null, null]});
    expect(totalsRow(big, tableColumns([fields[0]], []))).toEqual(['Total', '3,420']);
  });

  it('sums count columns only and leaves Project label and filter cells blank', () => {
    expect(totalsRow(projects, columns)).toEqual(['Total', '', '', '', '', '10', '14']);
    expect(totalsRow([projects[1]], columns)).toEqual(['Total', '', '', '', '', '', '']);
  });
});

describe('column visibility', () => {
  const columns = tableColumns(fields, []);

  it('lists toggleable columns without Project, in position order', () => {
    expect(labels(toggleableColumns(columns))).toEqual(['Species', 'Tissue', 'Completed', 'Funding Source', 'Subjects', 'Samples']);
  });

  it('hides chosen columns but never Project, keeping order', () => {
    expect(labels(visibleColumns(columns, ['tissue', 'subjects', 'project']))).toEqual(['Project', 'Species', 'Completed', 'Funding Source', 'Samples']);
    expect(visibleColumns(columns, [])).toEqual(columns);
  });

  it('parses stored ids, ignoring junk and ids that no longer exist', () => {
    expect(parseHiddenColumns('["tissue","gone","tissue"]', fields, [])).toEqual(['tissue']);
    expect(parseHiddenColumns('not json', fields, [])).toEqual([]);
    expect(parseHiddenColumns('{"a":1}', fields, [])).toEqual([]);
    expect(parseHiddenColumns(null, fields, [])).toEqual([]);
  });

  it('reads and writes the stored text', () => {
    const memory = (items) => ({
      items,
      getItem: (key) => (key in items ? items[key] : null),
      setItem: (key, value) => { items[key] = value; }
    });
    expect(readHiddenColumns(memory({[HIDDEN_COLUMNS_KEY]: '["species"]'}))).toEqual('["species"]');
    expect(readHiddenColumns(memory({}))).toBeNull();
    expect(readHiddenColumns({getItem: () => { throw new Error('denied'); }})).toBeNull();

    const store = memory({});
    writeHiddenColumns('["samples"]', store);
    expect(store.items).toEqual({[HIDDEN_COLUMNS_KEY]: '["samples"]'});
    expect(() => writeHiddenColumns('[]', {setItem: () => { throw new Error('full'); }})).not.toThrow();
  });

  it('toggles an id in and out', () => {
    expect(toggleHidden([], 'species')).toEqual(['species']);
    expect(toggleHidden(['species', 'tissue'], 'species')).toEqual(['tissue']);
  });

  it('applies to csv and totals', () => {
    const visible = visibleColumns(columns, ['species', 'samples']);
    expect(buildCsv(visible, csvRows(projects, visible)).split('\n')[0]).toEqual('Project,Tissue,Completed,Funding Source,Subjects');
    expect(totalsRow(projects, visible)).toEqual(['Total', '', '', '', '10']);
  });
});

describe('modality (derived from the model count)', () => {
  const modalities = [
    {id: 1, name: 'single cell', model: 'sc_seq', position: 1, query: ['sc_seq', '::count']},
    {id: 2, name: 'bulk', model: 'rna_seq', position: 2, query: ['rna_seq', '::count']},
    {id: 3, name: 'spatial', model: 'spatial', position: 3, query: ['spatial', '::count']}
  ];
  const afields = [
    {name: 'sp', label: 'Species', kind: 'filter', modalityId: null, query: []},
    {name: 'all', label: 'All subjects', kind: 'count', modalityId: null, query: []},
    {name: 'sc', label: 'Single cell pools', kind: 'count', modalityId: 1, query: []},
    {name: 'scs', label: 'Single cell subjects', kind: 'count', modalityId: 1, query: []},
    {name: 'bk', label: 'Bulk samples', kind: 'count', modalityId: 2, query: []}
  ];
  // field answers (5), then modality answers (3)
  const aresults = {
    a: ['human', 3, 2, 0, 0, 4, 0, null],
    b: ['rat', 5, 0, 0, null, 0, null, 0],
    c: ['mouse', 2, 0, 4, 1, 7, 2, null]
  };
  const aprojects = buildProjects(afields, modalities, aresults);

  it('derives modalities from the answers: a number > 0 only', () => {
    expect(aprojects.map((p) => p.values.modality)).toEqual([['single cell'], [], ['single cell', 'bulk']]);
    expect(projectModalities(modalities, [0, null, 'x'])).toEqual([]);
    expect(projectModalities(modalities, [1, 0, 3])).toEqual(['single cell', 'spatial']);
  });

  it('has no Modality facet or column without modalities', () => {
    const noModalities = buildProjects(afields, [], aresults);
    expect(noModalities[0].values.modality).toBeUndefined();
    expect(modalityFacet([])).toEqual([]);
    expect(tableColumns(afields, []).map((c) => c.label)).toEqual(['Project', 'Species', 'All subjects']);
  });

  it('lists every modality in position order with live counts, zero greyed', () => {
    expect(liveModalityOptions(modalities, aprojects, {}, afields)).toEqual([
      {value: 'single cell', count: 2, groups: []}, {value: 'bulk', count: 1, groups: []}, {value: 'spatial', count: 0, groups: []}
    ]);
    expect(liveModalityOptions(modalities, aprojects, {sp: ['human']}).map((o) => o.count)).toEqual([1, 0, 0]);
    expect(liveModalityOptions(modalities, aprojects, {modality: ['bulk']}).map((o) => o.count)).toEqual([2, 1, 0]);
  });

  it('filters by modality with OR', () => {
    expect(filterProjects(aprojects, {modality: ['single cell']}).map((p) => p.project_name)).toEqual(['a', 'c']);
    expect(filterProjects(aprojects, {modality: ['bulk', 'single cell']}).map((p) => p.project_name)).toEqual(['a', 'c']);
    expect(filterProjects(aprojects, {modality: ['single cell'], sp: ['rat']})).toEqual([]);
  });

  it('shows Modality after Project and only applicable counts', () => {
    expect(tableColumns(afields, modalities).map((c) => c.label)).toEqual(['Project', 'Modality', 'Species', 'All subjects']);
    expect(tableColumns(afields, modalities, ['bulk']).map((c) => c.label)).toEqual(['Project', 'Modality', 'Species', 'All subjects', 'Bulk samples']);
    expect(applicableCounts(afields, modalities, ['single cell', 'bulk']).map((f) => f.name)).toEqual(['all', 'sc', 'scs', 'bk']);
    const columns = tableColumns(afields, modalities, ['single cell']);
    expect(tableRows(aprojects, columns)[2]).toEqual(['c', 'single cell, bulk', 'mouse', '2', '0', '4']);
    expect(tableRows(aprojects, columns)[1][1]).toEqual('n/a');
    expect(totalsRow(aprojects, columns)).toEqual(['Total', '', '', '10', '2', '4']);
    expect(buildCsv(columns, csvRows(aprojects, columns)).split('\n')[0]).toEqual('Project,Modality,Species,All subjects,Single cell pools,Single cell subjects');
  });

  it('lets the Modality column be hidden', () => {
    expect(parseHiddenColumns('["modality","zzz"]', afields, modalities)).toEqual(['modality']);
    expect(parseHiddenColumns('["modality"]', afields, [])).toEqual([]);
    expect(visibleColumns(tableColumns(afields, modalities), ['modality']).map((c) => c.label)).toEqual(['Project', 'Species', 'All subjects']);
  });
});

describe('filters nested in the Modality facet', () => {
  const modalities = [{id: 1, name: 'single cell', position: 1, query: []}, {id: 2, name: 'bulk', position: 2, query: []}];
  const tfields = [
    {name: 'sp', label: 'Species', kind: 'filter', modalityId: null, query: []},
    {name: 'sub', label: 'Subtype', kind: 'filter', modalityId: 1, query: []},
    {name: 'sc', label: 'Pools', kind: 'count', modalityId: 1, query: []},
    {name: 'bk', label: 'Bulk', kind: 'count', modalityId: 2, query: []}
  ];
  // field answers (4), then modality answers (2)
  const tprojects = buildProjects(tfields, modalities, {
    a: ['human', ['RNA', 'CITE'], 2, 0, 1, 0],
    b: ['rat', ['ATAC'], 1, 0, 3, 0],
    c: ['mouse', null, 0, 3, 0, 9],
    d: ['human', ['RNA'], 4, 5, 2, 1]
  });
  const group = modalityGroup(tfields, modalities);
  const shownNames = (selected) => names(filterProjects(tprojects, selected, group));

  it('takes the nested filters out of the sidebar and into the Modality group', () => {
    expect(sidebarFacets(tfields, modalities).map((f) => f.name)).toEqual(['modality', 'sp']);
    expect(modalityFilterFields(tfields, modalities[0]).map((f) => f.name)).toEqual(['sub']);
    expect(modalityFilterFields(tfields, modalities[1])).toEqual([]);
    expect(group).toEqual(['modality', 'sub']);
    expect(modalityGroup(tfields, [])).toEqual([]);
  });

  it('always lists nested values with live counts, ignoring the Modality selection', () => {
    const expected = [
      {
        value: 'single cell',
        count: 3,
        groups: [{name: 'sub', label: 'Subtype', options: [
          {value: 'ATAC', count: 1}, {value: 'CITE', count: 1}, {value: 'RNA', count: 2}
        ]}]
      },
      {value: 'bulk', count: 2, groups: []}
    ];
    expect(liveModalityOptions(modalities, tprojects, {}, tfields)).toEqual(expected);
    expect(liveModalityOptions(modalities, tprojects, {sub: ['ATAC'], modality: ['bulk']}, tfields)).toEqual(expected);
    const human = liveModalityOptions(modalities, tprojects, {sp: ['human']}, tfields);
    expect(human.map((o) => o.count)).toEqual([2, 1]);
    expect(human[0].groups[0].options).toEqual([
      {value: 'ATAC', count: 0}, {value: 'CITE', count: 1}, {value: 'RNA', count: 2}
    ]);
  });

  it('ORs everything ticked in the Modality facet and ANDs it with other facets', () => {
    expect(shownNames({sub: ['ATAC']})).toEqual(['b']);
    expect(shownNames({sub: ['CITE', 'ATAC']})).toEqual(['a', 'b']);
    expect(shownNames({sub: ['ATAC'], modality: ['bulk']})).toEqual(['b', 'c', 'd']);
    expect(shownNames({modality: ['single cell']})).toEqual(['a', 'b', 'd']);
    expect(shownNames({sub: ['RNA'], sp: ['human']})).toEqual(['a', 'd']);
    expect(shownNames({sub: ['ATAC'], sp: ['human']})).toEqual([]);
    // without the group, facets AND as before
    expect(names(filterProjects(tprojects, {sub: ['ATAC'], modality: ['bulk']}))).toEqual([]);
  });

  it('counts other facets against the Modality group as one OR facet', () => {
    expect(liveFacetOptions('sp', tprojects, {sub: ['ATAC'], modality: ['bulk']}, group)).toEqual([
      {value: 'human', count: 1}, {value: 'mouse', count: 1}, {value: 'rat', count: 1}
    ]);
  });

  it('shows a modality\'s columns when it or one of its nested values is ticked', () => {
    expect(activeModalities({}, tfields, modalities)).toEqual([]);
    expect(activeModalities({sub: ['ATAC']}, tfields, modalities)).toEqual(['single cell']);
    expect(activeModalities({modality: ['bulk']}, tfields, modalities)).toEqual(['bulk']);
    expect(activeModalities({modality: ['bulk'], sub: ['RNA'], sp: ['rat']}, tfields, modalities)).toEqual(['single cell', 'bulk']);
    expect(labels(tableColumns(tfields, modalities, activeModalities({}, tfields, modalities)))).toEqual(['Project', 'Modality', 'Species']);
    expect(labels(tableColumns(tfields, modalities, activeModalities({sub: ['ATAC']}, tfields, modalities))))
      .toEqual(['Project', 'Modality', 'Species', 'Subtype', 'Pools']);
    expect(applicableFilters(tfields, modalities, ['bulk']).map((f) => f.name)).toEqual(['sp']);
  });

  it('shows the most specific modality: nested values, else the modality name', () => {
    expect(tprojects.map((p) => p.modalityDetail)).toEqual([['CITE', 'RNA'], ['ATAC'], ['bulk'], ['RNA', 'bulk']]);
    const columns = tableColumns(tfields, modalities, ['single cell']);
    expect(tableRows(tprojects, columns)[0]).toEqual(['a', 'CITE, RNA', 'human', 'CITE, RNA', '2']);
    expect(csvRows(tprojects, columns)[0][3]).toEqual('CITE; RNA');
    expect(tableRows(tprojects, columns)[2]).toEqual(['c', 'bulk', 'mouse', 'n/a', '0']);
    expect(csvRows(tprojects, columns)[3][1]).toEqual('RNA; bulk');
    expect(modalityDetail({values: {modality: [], sub: null}}, tfields, modalities)).toEqual([]);
    expect(tableRows([{project_name: 'e', values: {modality: []}, modalityDetail: []}], columns.slice(0, 2))).toEqual([['e', 'n/a']]);
  });
});
