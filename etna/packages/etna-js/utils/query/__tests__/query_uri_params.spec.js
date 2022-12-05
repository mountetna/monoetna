"use strict";
exports.__esModule = true;
var query_uri_params_1 = require("../query_uri_params");
describe('migrateSubclauses', function () {
    describe('old clause format', function () {
        it('correctly migrates', function () {
            var input = [
                {
                    anyMap: {},
                    modelName: 'labor',
                    clauses: [
                        {
                            modelName: 'labor',
                            any: true,
                            attributeName: 'name',
                            attributeType: 'string',
                            operator: '::equals',
                            operand: 'foo'
                        }
                    ]
                }
            ];
            expect(query_uri_params_1.migrateSubclauses(input)).toEqual([
                {
                    anyMap: {},
                    modelName: 'labor',
                    clauses: [
                        {
                            modelName: 'labor',
                            any: true,
                            subclauses: [
                                {
                                    attributeName: 'name',
                                    attributeType: 'string',
                                    operator: '::equals',
                                    operand: 'foo'
                                }
                            ]
                        }
                    ]
                }
            ]);
        });
        it('correctly migrates with no operand', function () {
            var input = [
                {
                    anyMap: {},
                    modelName: 'labor',
                    clauses: [
                        {
                            modelName: 'labor',
                            any: true,
                            attributeName: 'name',
                            attributeType: 'string',
                            operator: '::has',
                            operand: ''
                        }
                    ]
                }
            ];
            expect(query_uri_params_1.migrateSubclauses(input)).toEqual([
                {
                    anyMap: {},
                    modelName: 'labor',
                    clauses: [
                        {
                            modelName: 'labor',
                            any: true,
                            subclauses: [
                                {
                                    attributeName: 'name',
                                    attributeType: 'string',
                                    operator: '::has',
                                    operand: ''
                                }
                            ]
                        }
                    ]
                }
            ]);
        });
    });
    it('generates empty subclause with empty old format', function () {
        var input = [
            {
                anyMap: {},
                modelName: 'labor',
                clauses: [
                    {
                        modelName: 'labor',
                        any: true,
                        attributeName: '',
                        attributeType: '',
                        operator: '',
                        operand: ''
                    }
                ]
            }
        ];
        expect(query_uri_params_1.migrateSubclauses(input)).toEqual([
            {
                anyMap: {},
                modelName: 'labor',
                clauses: [
                    {
                        modelName: 'labor',
                        any: true,
                        subclauses: [
                            {
                                attributeName: '',
                                attributeType: '',
                                operator: '',
                                operand: ''
                            }
                        ]
                    }
                ]
            }
        ]);
    });
    it('ignores new clause format', function () {
        var input = [
            {
                anyMap: {},
                modelName: 'labor',
                clauses: [
                    {
                        modelName: 'labor',
                        any: true,
                        subclauses: [
                            {
                                attributeName: 'name',
                                attributeType: 'string',
                                operator: '::equals',
                                operand: 'foo'
                            }
                        ]
                    }
                ]
            }
        ];
        expect(query_uri_params_1.migrateSubclauses(input)).toEqual([
            {
                anyMap: {},
                modelName: 'labor',
                clauses: [
                    {
                        modelName: 'labor',
                        any: true,
                        subclauses: [
                            {
                                attributeName: 'name',
                                attributeType: 'string',
                                operator: '::equals',
                                operand: 'foo'
                            }
                        ]
                    }
                ]
            }
        ]);
    });
});
