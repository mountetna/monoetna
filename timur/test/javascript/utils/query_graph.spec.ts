import {QueryGraph} from '../../../lib/client/jsx/utils/query_graph';
import {models} from '../fixtures/models';

describe('QueryGraph', () => {
  let graph: QueryGraph;

  beforeEach(() => {
    graph = new QueryGraph(models);
  });

  it('includes the project model', () => {
    expect(Object.keys(graph.graph.children).includes('project')).toEqual(true);
    expect(Object.keys(graph.graph.parents).includes('project')).toEqual(true);
  });

  it('adds table and link connections', () => {
    expect(Object.keys(graph.graph.children).includes('prize')).toEqual(true);
    expect(Object.keys(graph.graph.parents).includes('prize')).toEqual(true);
    expect(graph.pathsFrom('labor')).toEqual([
      ['labor', 'monster', 'victim', 'demographics'],
      ['labor', 'monster', 'habitat', 'vegetation'],
      ['labor', 'monster', 'victim', 'wound'],
      ['labor', 'monster', 'aspect'],
      ['labor', 'prize']
    ]);
  });

  it('provides all paths from a child model, up and down the graph', () => {
    expect(graph.allPaths('prize')).toEqual([
      ['labor', 'project'],
      ['labor', 'monster', 'victim', 'demographics'],
      ['labor', 'monster', 'habitat', 'vegetation'],
      ['labor', 'monster', 'victim', 'wound'],
      ['labor', 'monster', 'aspect']
    ]);
  });

  it('correctly returns children and one-to-many status', () => {
    expect(graph.childrenMap('monster')).toEqual({
      aspect: true,
      habitat: false,
      victim: true
    });

    expect(graph.childrenMap('habitat')).toEqual({
      vegetation: true
    });

    expect(graph.childrenMap('wound')).toEqual({});
  });

  describe('sliceable', () => {
    it('does not allow root models to be sliced', async () => {
      expect(graph.sliceable('prize','prize')).toEqual(false);
    });

    it('allows non-root tables to be sliced', async () => {
      expect(graph.sliceable('monster','monster')).toEqual(false);
      expect(graph.sliceable('monster','prize')).toEqual(true);
      expect(graph.sliceable('monster','aspect')).toEqual(true);
    });
  });

  describe('for xcrs1 models', () => {
    const models = require('../fixtures/xcrs1_magma_metadata.json').models;
    beforeEach(() => {
      graph = new QueryGraph(models);
    });

    it('handles the path laterally from subject -> sc_seq', () => {
      expect(graph.shortestPath('subject', 'sc_seq')).toEqual([
        'experiment',
        'biospecimen_group',
        'sc_seq'
      ]);
    });

    it('handles the path laterally from sc_seq -> subject', () => {
      expect(graph.shortestPath('sc_seq', 'subject')).toEqual([
        'biospecimen_group',
	'experiment',
        'subject'
      ]);
    });

    it('handles the path laterally from sc_seq -> cytof', () => {
      expect(graph.shortestPath('sc_seq', 'cytof')).toEqual([
        'biospecimen_group',
        'biospecimen',
        'cytof'
      ]);
    });
  });
});
