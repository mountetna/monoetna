"use strict";
exports.__esModule = true;
var directed_graph_1 = require("../directed_graph");
describe('DirectedGraph', function () {
    var graphOne;
    var graphTwo;
    var graphThree;
    beforeEach(function () {
        graphOne = new directed_graph_1.DirectedGraph();
        graphOne.addConnection('a', 'b');
        graphOne.addConnection('a', 'c');
        graphOne.addConnection('b', 'd');
        graphOne.addConnection('b', 'e');
        graphOne.addConnection('c', 'f');
        graphOne.addConnection('e', 'g');
        graphTwo = new directed_graph_1.DirectedGraph();
        graphTwo.addConnection('a', 'b');
        graphTwo.addConnection('a', 'c');
        graphTwo.addConnection('c', 'd');
        graphTwo.addConnection('b', 'd');
        graphTwo.addConnection('d', 'e');
        graphTwo.addConnection('d', 'f');
        graphTwo.addConnection('e', 'g');
        graphTwo.addConnection('f', 'g');
        graphTwo.addConnection('g', 'h');
        graphTwo.addConnection('g', 'i');
        graphTwo.addConnection('t', 'i');
        graphThree = new directed_graph_1.DirectedGraph();
        graphThree.addConnection('a', 'b');
        graphThree.addConnection('a', 'c');
        graphThree.addConnection('c', 'd');
        graphThree.addConnection('d', 'b');
        graphThree.addConnection('b', 'e');
        graphThree.addConnection('a', 'f');
        graphThree.addConnection('e', 'f');
    });
    describe('pathsFrom', function () {
        it('divergent paths work', function () {
            expect(graphOne.pathsFrom('a')).toEqual([
                ['a', 'b', 'e', 'g'],
                ['a', 'b', 'd'],
                ['a', 'c', 'f']
            ]);
            expect(graphOne.pathsFrom('a', false)).toEqual([
                ['b', 'e', 'g'],
                ['b', 'd'],
                ['c', 'f']
            ]);
        });
        it('convergent and divergent paths work', function () {
            expect(graphTwo.pathsFrom('a')).toEqual([
                ['a', 'c', 'd', 'f', 'g', 'h'],
                ['a', 'c', 'd', 'f', 'g', 'i'],
                ['a', 'c', 'd', 'e', 'g'],
                ['a', 'b', 'd']
            ]);
        });
    });
    describe('#serializedPathFrom', function () {
        it('generates a linear serialization of convergent and divergent paths', function () {
            expect(graphTwo.serializedPathFrom('a')).toEqual([
                'a',
                'c',
                'b',
                'd',
                'f',
                'e',
                'g',
                'h',
                'i'
            ]);
            expect(graphTwo.serializedPathFrom('a', false)).toEqual([
                'c',
                'b',
                'd',
                'f',
                'e',
                'g',
                'h',
                'i'
            ]);
        });
    });
    describe('#asNormalizedHash', function () {
        it('works', function () {
            expect(graphTwo.asNormalizedHash('a')).toEqual({
                a: ['b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'],
                b: ['d', 'e', 'f', 'g', 'h', 'i'],
                c: ['d', 'e', 'f', 'g', 'h', 'i'],
                d: ['e', 'f', 'g', 'h', 'i'],
                e: ['g', 'h', 'i'],
                f: ['g', 'h', 'i'],
                g: ['h', 'i'],
                h: [],
                i: []
            });
            expect(graphTwo.asNormalizedHash('a', false)).toEqual({
                b: ['d', 'e', 'f', 'g', 'h', 'i'],
                c: ['d', 'e', 'f', 'g', 'h', 'i'],
                d: ['e', 'f', 'g', 'h', 'i'],
                e: ['g', 'h', 'i'],
                f: ['g', 'h', 'i'],
                g: ['h', 'i'],
                h: [],
                i: []
            });
            expect(graphThree.asNormalizedHash('a')).toEqual({
                a: ['b', 'c', 'f', 'e', 'd'],
                b: ['e', 'f'],
                c: ['e', 'd', 'f', 'b'],
                d: ['f', 'b', 'e'],
                e: ['f'],
                f: []
            });
        });
    });
});
