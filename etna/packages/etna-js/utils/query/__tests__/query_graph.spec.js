"use strict";
exports.__esModule = true;
var query_graph_1 = require("../query_graph");
var models = {
    monster: {
        documents: {},
        revisions: {},
        views: {},
        template: require('../fixtures/template_monster.json')
    },
    labor: {
        documents: {},
        revisions: {},
        views: {},
        template: require('../fixtures/template_labor.json')
    },
    project: {
        documents: {},
        revisions: {},
        views: {},
        template: require('../fixtures/template_project.json')
    },
    victim: {
        documents: {},
        revisions: {},
        views: {},
        template: require('../fixtures/template_victim.json')
    },
    wound: {
        documents: {},
        revisions: {},
        views: {},
        template: require('../fixtures/template_wound.json')
    },
    habitat: {
        documents: {},
        revisions: {},
        views: {},
        template: require('../fixtures/template_habitat.json')
    },
    vegetation: {
        documents: {},
        revisions: {},
        views: {},
        template: require('../fixtures/template_vegetation.json')
    }
};
describe('QueryGraph', function () {
    var graph;
    beforeEach(function () {
        graph = new query_graph_1.QueryGraph(models);
    });
    it('ignores the project model', function () {
        expect(Object.keys(graph.graph.children).includes('project')).toEqual(false);
        expect(Object.keys(graph.graph.parents).includes('project')).toEqual(false);
    });
    it('adds table and link connections', function () {
        expect(Object.keys(graph.graph.children).includes('prize')).toEqual(true);
        expect(Object.keys(graph.graph.parents).includes('prize')).toEqual(true);
        expect(graph.pathsFrom('labor')).toEqual([
            ['labor', 'monster', 'habitat', 'vegetation'],
            ['labor', 'monster', 'victim', 'wound'],
            ['labor', 'prize']
        ]);
    });
    it('provides all paths from a child model, up and down the graph', function () {
        expect(graph.allPaths('prize')).toEqual([
            ['labor'],
            ['labor', 'monster', 'habitat', 'vegetation'],
            ['labor', 'monster', 'victim', 'wound'],
            ['labor', 'prize']
        ]);
    });
    it('correctly returns children and one-to-many status', function () {
        expect(graph.childrenMap('monster')).toEqual({
            habitat: false,
            victim: true,
            monster: false
        });
        expect(graph.childrenMap('habitat')).toEqual({
            vegetation: true,
            habitat: false
        });
        expect(graph.childrenMap('wound')).toEqual({
            wound: false
        });
    });
    describe('for xcrs1 models', function () {
        var models = require('../fixtures/xcrs1_magma_metadata.json').models;
        beforeEach(function () {
            graph = new query_graph_1.QueryGraph(models);
        });
        it('handles the path laterally from subject -> sc_seq', function () {
            expect(graph.shortestPath('subject', 'sc_seq')).toEqual([
                'biospecimen',
                'biospecimen_group',
                'sc_seq'
            ]);
        });
        it('handles the path laterally from sc_seq -> subject', function () {
            expect(graph.shortestPath('sc_seq', 'subject')).toEqual([
                'biospecimen_group',
                'biospecimen',
                'subject'
            ]);
        });
        it('handles the path laterally from sc_seq -> cytof', function () {
            expect(graph.shortestPath('sc_seq', 'cytof')).toEqual([
                'biospecimen_group',
                'biospecimen',
                'cytof'
            ]);
        });
    });
});
