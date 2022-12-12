"use strict";
exports.__esModule = true;
var query_deserializer_1 = require("../query_deserializer");
describe('QueryDeserializer', function () {
    it('identifies the raw filters', function () {
        var deserializer = new query_deserializer_1.QueryDeserializer([
            'monster',
            [
                '::and',
                ['labor', ['name', '::in', ['lion', 'hydra', 'apples']], '::any'],
                ['name', '::equals', 'Nemean Lion'],
                ['labor', ['number', '::equals', 2], '::any'],
                ['labor', ['prize', ['name', '::equals', 'Apples'], '::any'], '::any']
            ],
            '::all',
            [
                ['species'],
                ['stats', '::url'],
                ['labor', 'year'],
                ['labor', 'completed'],
                ['labor', 'contributions', '::slice', ['Athens', 'Sidon']],
                [
                    'labor',
                    'prize',
                    ['name', '::equals', 'Sparta'],
                    '::first',
                    'value'
                ],
                ['victim', '::first', 'country']
            ]
        ], [
            'monster.species',
            'monster.stats',
            'labor.year',
            'labor.constributions.Athens',
            'labor.contributions.Sidon',
            'prize.value',
            'victim.country'
        ]);
        expect(deserializer.rawFilters()).toEqual([
            '::and',
            ['labor', ['name', '::in', ['lion', 'hydra', 'apples']], '::any'],
            ['name', '::equals', 'Nemean Lion'],
            ['labor', ['number', '::equals', 2], '::any'],
            ['labor', ['prize', ['name', '::equals', 'Apples'], '::any'], '::any']
        ]);
    });
    it('identifies the raw columns', function () {
        var deserializer = new query_deserializer_1.QueryDeserializer([
            'monster',
            [
                '::and',
                ['labor', ['name', '::in', ['lion', 'hydra', 'apples']], '::any'],
                ['name', '::equals', 'Nemean Lion'],
                ['labor', ['number', '::equals', 2], '::any'],
                ['labor', ['prize', ['name', '::equals', 'Apples'], '::any'], '::any']
            ],
            '::all',
            [
                ['species'],
                ['stats', '::url'],
                ['labor', 'year'],
                ['labor', 'completed'],
                ['labor', 'contributions', '::slice', ['Athens', 'Sidon']],
                [
                    'labor',
                    'prize',
                    ['name', '::equals', 'Sparta'],
                    '::first',
                    'value'
                ],
                ['victim', '::first', 'country']
            ]
        ], [
            'monster.species',
            'monster.stats',
            'labor.year',
            'labor.constributions.Athens',
            'labor.contributions.Sidon',
            'prize.value',
            'victim.country'
        ]);
        expect(deserializer.rawColumns()).toEqual([
            ['species'],
            ['stats', '::url'],
            ['labor', 'year'],
            ['labor', 'completed'],
            ['labor', 'contributions', '::slice', ['Athens', 'Sidon']],
            ['labor', 'prize', ['name', '::equals', 'Sparta'], '::first', 'value'],
            ['victim', '::first', 'country']
        ]);
    });
    describe('orRecordFilterIndices', function () {
        it('gives them correct indices when no ::and filters', function () {
            var deserializer = new query_deserializer_1.QueryDeserializer([
                'monster',
                [
                    '::and',
                    [
                        '::or',
                        ['labor', ['name', '::in', ['lion', 'hydra', 'apples']], '::any'],
                        ['labor', ['number', '::equals', 2], '::any']
                    ]
                ],
                '::all',
                [
                    ['species'],
                    ['stats', '::url'],
                    ['labor', 'year'],
                    ['labor', 'completed'],
                    ['labor', 'contributions', '::slice', ['Athens', 'Sidon']],
                    [
                        'labor',
                        'prize',
                        ['name', '::equals', 'Sparta'],
                        '::all',
                        'value'
                    ],
                    ['victim', '::all', 'country']
                ]
            ], [
                'monster.species',
                'monster.stats',
                'labor.year',
                'labor.constributions.Athens',
                'labor.contributions.Sidon',
                'prize.value',
                'victim.country'
            ]);
            expect(deserializer.orRecordFilterIndices()).toEqual([0, 1]);
        });
        it('gives them correct indices when at the end of filter list', function () {
            var deserializer = new query_deserializer_1.QueryDeserializer([
                'monster',
                [
                    '::and',
                    ['name', '::equals', 'Nemean Lion'],
                    [
                        'labor',
                        ['prize', ['name', '::equals', 'Apples'], '::any'],
                        '::any'
                    ],
                    [
                        '::or',
                        ['labor', ['name', '::in', ['lion', 'hydra', 'apples']], '::any'],
                        ['labor', ['number', '::equals', 2], '::any']
                    ]
                ],
                '::all',
                [
                    ['species'],
                    ['stats', '::url'],
                    ['labor', 'year'],
                    ['labor', 'completed'],
                    ['labor', 'contributions', '::slice', ['Athens', 'Sidon']],
                    [
                        'labor',
                        'prize',
                        ['name', '::equals', 'Sparta'],
                        '::all',
                        'value'
                    ],
                    ['victim', '::all', 'country']
                ]
            ], [
                'monster.species',
                'monster.stats',
                'labor.year',
                'labor.constributions.Athens',
                'labor.contributions.Sidon',
                'prize.value',
                'victim.country'
            ]);
            expect(deserializer.orRecordFilterIndices()).toEqual([2, 3]);
        });
    });
    it('identifies rootModel', function () {
        var deserializer = new query_deserializer_1.QueryDeserializer([
            'monster',
            [
                '::and',
                ['labor', ['name', '::in', ['lion', 'hydra', 'apples']], '::any'],
                ['name', '::equals', 'Nemean Lion']
            ],
            '::all',
            ['name']
        ], ['monster.name']);
        expect(deserializer.rootModel()).toEqual('monster');
    });
    describe('#flattenFilter', function () {
        fit('handles nested ::and and ::or filters', function () {
            var deserializer = new query_deserializer_1.QueryDeserializer([
                'monster',
                [
                    '::and',
                    ['labor', ['name', '::in', ['lion', 'hydra', 'apples']], '::any'],
                    ['name', '::equals', 'Nemean Lion'],
                    [
                        '::or',
                        ['labor', ['year', '::<', 2022], '::every'],
                        ['labor', ['number', '::equals', 2], '::any']
                    ]
                ],
                '::all',
                ['name']
            ], ['monster.name']);
            expect(deserializer.flattenFilters([
                '::and',
                ['labor', ['name', '::in', ['lion', 'hydra', 'apples']], '::any'],
                ['name', '::equals', 'Nemean Lion'],
                [
                    '::or',
                    ['labor', ['year', '::<', 2022], '::every'],
                    ['labor', ['number', '::equals', 2], '::any']
                ]
            ])).toEqual([
                {
                    anyMap: {},
                    clauses: [
                        {
                            any: true,
                            modelName: 'labor',
                            subclauses: [
                                {
                                    attributeName: 'name',
                                    attributeType: '',
                                    operand: 'lion,hydra,apples',
                                    operator: '::in'
                                }
                            ]
                        }
                    ],
                    modelName: 'labor'
                },
                {
                    anyMap: {},
                    clauses: [
                        {
                            any: true,
                            modelName: 'monster',
                            subclauses: [
                                {
                                    attributeName: 'name',
                                    attributeType: '',
                                    operand: 'Nemean Lion',
                                    operator: '::equals'
                                }
                            ]
                        }
                    ],
                    modelName: 'monster'
                },
                {
                    anyMap: {
                        labor: false
                    },
                    clauses: [
                        {
                            any: true,
                            modelName: 'labor',
                            subclauses: [
                                {
                                    attributeName: 'year',
                                    attributeType: '',
                                    operand: 2022,
                                    operator: '::<'
                                }
                            ]
                        }
                    ],
                    modelName: 'labor'
                },
                {
                    anyMap: {
                        labor: false
                    },
                    clauses: [
                        {
                            any: true,
                            modelName: 'labor',
                            subclauses: [
                                {
                                    attributeName: 'number',
                                    attributeType: '',
                                    operand: 2,
                                    operator: '::equals'
                                }
                            ]
                        }
                    ],
                    modelName: 'labor'
                }
            ]);
        });
    });
    describe('recordFilters', function () {
        it('works with just one on root model', function () {
            var deserializer = new query_deserializer_1.QueryDeserializer(['monster', ['name', '::equals', 'Nemean Lion'], '::all', ['name']], ['monster.name']);
            expect(deserializer.recordFilters()).toEqual([
                {
                    modelName: 'monster',
                    anyMap: {},
                    clauses: [
                        {
                            subclauses: [
                                {
                                    attributeName: 'name',
                                    operator: '::equals',
                                    operand: 'Nemean Lion',
                                    attributeType: '' // Hopefully this gets populated correctly on the UI side, later...
                                }
                            ],
                            modelName: 'monster',
                            any: true
                        }
                    ]
                }
            ]);
        });
        it('works with multiple model filters joined by ::and', function () { });
        it('correctly sets ::every in anyMap', function () { });
        it('includes multiple subclauses in root model filter', function () { });
        it('includes multiple subclauses in non-root model filter', function () { });
        it('comma-joins ::in operands', function () {
            var deserializer = new query_deserializer_1.QueryDeserializer([
                'monster',
                ['name', '::in', ['Nemean Lion', 'Lernean Hydra']],
                '::all',
                ['name']
            ], ['monster.name']);
            expect(deserializer.recordFilters()).toEqual([
                {
                    modelName: 'monster',
                    anyMap: {},
                    clauses: [
                        {
                            subclauses: [
                                {
                                    attributeName: 'name',
                                    operator: '::in',
                                    operand: 'Nemean Lion,Lernean Hydra',
                                    attributeType: '' // Hopefully this gets populated correctly on the UI side, later...
                                }
                            ],
                            modelName: 'monster',
                            any: true
                        }
                    ]
                }
            ]);
        });
        it('correctly constructs inverted operators', function () { });
        it('correctly constructs terminal operators', function () { });
    });
    describe('columns', function () {
        it('works with no edited user_columns', function () { });
        it('works with edited user_columns', function () { });
        it('comma-joins ::slice operands', function () { });
    });
});
