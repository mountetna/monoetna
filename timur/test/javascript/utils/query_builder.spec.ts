import {QueryBuilder} from '../../../lib/client/jsx/utils/query_builder';
import {QueryGraph} from '../../../lib/client/jsx/utils/query_graph';
import {QuerySlice} from '../../../lib/client/jsx/contexts/query/query_types';

function stampTemplate(template: any) {
  return {
    documents: {},
    revisions: {},
    views: {},
    template
  };
}

const models = {
  monster: stampTemplate(require('../fixtures/template_monster.json')),
  labor: stampTemplate(require('../fixtures/template_labor.json')),
  project: stampTemplate(require('../fixtures/template_project.json')),
  prize: stampTemplate(require('../fixtures/template_prize.json')),
  victim: stampTemplate(require('../fixtures/template_victim.json')),
  wound: stampTemplate(require('../fixtures/template_wound.json'))
};

describe('QueryBuilder', () => {
  let graph: QueryGraph;
  let builder: QueryBuilder;

  beforeEach(() => {
    graph = new QueryGraph(models);
    builder = new QueryBuilder(graph);
  });

  function stamp(
    model_name: string,
    attribute_name: string,
    slices: QuerySlice[]
  ) {
    return {
      model_name,
      attribute_name,
      display_label: `${model_name}.${attribute_name}`,
      slices: slices
    };
  }

  describe('for xcrs1 models', () => {
    const models = require('../fixtures/xcrs1_magma_metadata.json').models;
    beforeEach(() => {
      graph = new QueryGraph(models);
      builder = new QueryBuilder(graph);
    });

    it('handles biospecimen -> sc_seq', () => {
      builder.addRootModel('subject');
      builder.addRecordFilters([
        {
          anyMap: {biospecimen: true, sc_seq: true},
          clauses: [
            {
              subclauses: [
                {
                  attributeName: 'tube_name',
                  operand: '',
                  operator: '::has',
                  attributeType: 'text'
                }
              ],
              modelName: 'sc_seq',
              any: true
            }
          ],
          modelName: 'sc_seq'
        }
      ]);
      builder.addColumns([stamp('subject', 'name', [])]);
      expect(builder.query()).toEqual([
        'subject',
        [
          'experiment',
          [
            'biospecimen_group',
            ['sc_seq', ['::has', 'tube_name'], '::any'],
            '::any'
          ],
          '::any'
        ],
        '::all',
        ['name']
      ]);
    });
  });

  it('works', () => {
    builder.addRootModel('monster');
    builder.addColumns([
      stamp('monster', 'name', []),
      stamp('monster', 'species', []),
      stamp('monster', 'stats', []),
      stamp('labor', 'year', []),
      stamp('labor', 'completed', []),
      stamp('labor', 'contributions', [
        {
          modelName: 'labor',
          clause: {
            subclauses: [
              {
                attributeName: 'contributions',
                operator: '::slice',
                operand: 'Athens,Sidon',
                attributeType: 'matrix'
              }
            ],
            modelName: 'labor',
            any: true
          }
        }
      ]),
      stamp('prize', 'value', [
        {
          modelName: 'prize',
          clause: {
            subclauses: [
              {
                attributeName: 'name',
                operator: '::equals',
                operand: 'Sparta',
                attributeType: 'text'
              }
            ],
            modelName: 'prize',
            any: true
          }
        }
      ]),
      stamp('victim', 'country', [])
    ]);
    builder.addRecordFilters([
      {
        modelName: 'labor',
        anyMap: {},
        clauses: [
          {
            subclauses: [
              {
                attributeName: 'name',
                operator: '::in',
                operand: 'lion,hydra,apples',
                attributeType: 'text'
              }
            ],
            modelName: 'labor',
            any: true
          }
        ]
      },
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
                attributeType: 'text'
              }
            ],
            modelName: 'monster',
            any: true
          }
        ]
      },
      {
        modelName: 'labor',
        anyMap: {},
        clauses: [
          {
            subclauses: [
              {
                attributeName: 'number',
                operator: '::equals',
                operand: '2',
                attributeType: 'number'
              }
            ],
            modelName: 'labor',
            any: true
          }
        ]
      },
      {
        modelName: 'prize',
        clauses: [
          {
            subclauses: [
              {
                attributeName: 'name',
                operator: '::equals',
                operand: 'Apples',
                attributeType: 'text'
              }
            ],
            modelName: 'prize',
            any: true
          }
        ],
        anyMap: {
          prize: true
        }
      }
    ]);

    expect(builder.query()).toEqual([
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
        ['labor', 'prize', ['name', '::equals', 'Sparta'], '::first', 'value'],
        ['victim', '::first', 'country']
      ]
    ]);

    builder.setFlatten(false);

    expect(builder.query()).toEqual([
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
        ['labor', 'prize', ['name', '::equals', 'Sparta'], '::all', 'value'],
        ['victim', '::all', 'country']
      ]
    ]);

    builder.setOrRecordFilterIndices([0, 2]);

    expect(builder.query()).toEqual([
      'monster',
      [
        '::and',
        ['name', '::equals', 'Nemean Lion'],
        ['labor', ['prize', ['name', '::equals', 'Apples'], '::any'], '::any'],
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
        ['labor', 'prize', ['name', '::equals', 'Sparta'], '::all', 'value'],
        ['victim', '::all', 'country']
      ]
    ]);
  });

  it('combines multiple clauses in filters', () => {
    builder.addRootModel('monster');
    builder.addColumns([stamp('monster', 'name', [])]);
    builder.addRecordFilters([
      {
        modelName: 'labor',
        anyMap: {},
        clauses: [
          {
            subclauses: [
              {
                attributeName: 'name',
                operator: '::in',
                operand: 'lion,hydra,apples',
                attributeType: 'text'
              }
            ],
            modelName: 'labor',
            any: true
          },
          {
            subclauses: [
              {
                attributeName: 'number',
                operator: '::>',
                operand: '2',
                attributeType: 'number'
              }
            ],
            modelName: 'labor',
            any: true
          },
          {
            subclauses: [
              {
                attributeName: 'number',
                operator: '::<=',
                operand: '8',
                attributeType: 'number'
              }
            ],
            modelName: 'labor',
            any: true
          }
        ]
      }
    ]);

    expect(builder.query()).toEqual([
      'monster',
      [
        'labor',
        [
          '::and',
          ['name', '::in', ['lion', 'hydra', 'apples']],
          ['number', '::>', 2],
          ['number', '::<=', 8]
        ],
        '::any'
      ],
      '::all',
      ['name']
    ]);
  });

  it('combines multiple subclauses in filters', () => {
    builder.addRootModel('monster');
    builder.addColumns([stamp('monster', 'name', [])]);
    builder.addRecordFilters([
      {
        modelName: 'labor',
        anyMap: {},
        clauses: [
          {
            subclauses: [
              {
                attributeName: 'name',
                operator: '::in',
                operand: 'lion,hydra,apples',
                attributeType: 'text'
              }
            ],
            modelName: 'labor',
            any: true
          },
          {
            subclauses: [
              {
                attributeName: 'number',
                operator: '::>',
                operand: '2',
                attributeType: 'number'
              },
              {
                attributeName: 'number',
                operator: '::<=',
                operand: '8',
                attributeType: 'number'
              }
            ],
            modelName: 'labor',
            any: true
          }
        ]
      }
    ]);

    expect(builder.query()).toEqual([
      'monster',
      [
        'labor',
        [
          '::and',
          ['name', '::in', ['lion', 'hydra', 'apples']],
          ['::and', ['number', '::>', 2], ['number', '::<=', 8]]
        ],
        '::any'
      ],
      '::all',
      ['name']
    ]);
  });

  it('does not wrap individual filters on root model in ::and', () => {
    builder.addRootModel('monster');
    builder.addColumns([stamp('monster', 'name', [])]);
    builder.addRecordFilters([
      {
        modelName: 'monster',
        anyMap: {},
        clauses: [
          {
            subclauses: [
              {
                attributeName: 'name',
                operator: '::in',
                operand: 'lion,hydra,apples',
                attributeType: 'text'
              }
            ],
            modelName: 'monster',
            any: true
          }
        ]
      }
    ]);

    expect(builder.query()).toEqual([
      'monster',
      ['name', '::in', ['lion', 'hydra', 'apples']],
      '::all',
      ['name']
    ]);
  });

  it('adds slice for root model', () => {
    builder.addRootModel('labor');
    builder.addColumns([
      stamp('labor', 'name', []),
      stamp('labor', 'contributions', [
        {
          modelName: 'labor',
          clause: {
            subclauses: [
              {
                attributeName: 'contributions',
                operator: '::slice',
                operand: 'Athens,Sidon',
                attributeType: 'matrix'
              }
            ],
            modelName: 'labor',
            any: true
          }
        }
      ])
    ]);

    expect(builder.query()).toEqual([
      'labor',
      '::all',
      [['contributions', '::slice', ['Athens', 'Sidon']]]
    ]);
  });

  it('adds matrix slice with multiple clauses', () => {
    builder.addRootModel('labor');
    builder.addColumns([
      stamp('labor', 'name', []),
      stamp('victim', 'weapons', [
        {
          modelName: 'victim',
          clause: {
            subclauses: [
              {
                attributeName: 'weapons',
                operator: '::slice',
                operand: 'sword,hands',
                attributeType: 'matrix'
              }
            ],
            modelName: 'victim',
            any: true
          }
        },
        {
          modelName: 'monster',
          clause: {
            subclauses: [
              {
                attributeName: 'name',
                operator: '::matches',
                operand: 'ion',
                attributeType: 'text'
              }
            ],
            modelName: 'monster',
            any: true
          }
        }
      ])
    ]);

    expect(builder.query()).toEqual([
      'labor',
      '::all',
      [
        [
          'monster',
          ['name', '::matches', 'ion'],
          'victim',
          '::first',
          'weapons',
          '::slice',
          ['sword', 'hands']
        ]
      ]
    ]);
  });

  it('returns a count query string', () => {
    builder.addRootModel('monster');
    builder.addColumns([
      stamp('monster', 'name', []),
      stamp('labor', 'year', []),
      stamp('labor', 'completed', []),
      stamp('monster', 'species', []),
      stamp('prize', 'value', [
        {
          modelName: 'prize',
          clause: {
            subclauses: [
              {
                attributeName: 'name',
                operator: '::equals',
                operand: 'Sparta',
                attributeType: 'text'
              }
            ],
            modelName: 'prize',
            any: true
          }
        }
      ])
    ]);
    builder.addRecordFilters([
      {
        modelName: 'labor',
        clauses: [
          {
            subclauses: [
              {
                attributeName: 'name',
                operator: '::in',
                operand: 'lion,hydra,apples',
                attributeType: 'text'
              }
            ],
            modelName: 'labor',
            any: true
          }
        ],
        anyMap: {}
      },
      {
        modelName: 'monster',
        clauses: [
          {
            subclauses: [
              {
                attributeName: 'name',
                operator: '::equals',
                operand: 'Nemean Lion',
                attributeType: 'text'
              }
            ],
            modelName: 'monster',
            any: true
          }
        ],
        anyMap: {}
      }
    ]);

    expect(builder.count()).toEqual([
      'monster',
      [
        '::and',
        ['labor', ['name', '::in', ['lion', 'hydra', 'apples']], '::any'],
        ['name', '::equals', 'Nemean Lion']
      ],
      '::count'
    ]);
  });

  describe('handles any for', () => {
    it('deep paths in filters with some non-branching models', () => {
      builder.addRootModel('labor');
      builder.addRecordFilters([
        {
          modelName: 'wound',
          clauses: [
            {
              subclauses: [
                {
                  attributeName: 'location',
                  operator: '::equals',
                  operand: 'arm',
                  attributeType: 'text'
                }
              ],
              modelName: 'wound',
              any: true
            }
          ],
          anyMap: {
            victim: true,
            wound: true
          }
        }
      ]);
      builder.addColumns([stamp('labor', 'name', [])]);

      expect(builder.query()).toEqual([
        'labor',
        [
          'monster',
          [
            'victim',
            ['wound', ['location', '::equals', 'arm'], '::any'],
            '::any'
          ],
          '::any'
        ],
        '::all',
        ['name']
      ]);
    });

    it('deep paths in filters with branching models only', () => {
      builder.addRootModel('monster');
      builder.addRecordFilters([
        {
          modelName: 'wound',
          clauses: [
            {
              subclauses: [
                {
                  attributeName: 'location',
                  operator: '::equals',
                  operand: 'arm',
                  attributeType: 'text'
                }
              ],
              modelName: 'wound',
              any: true
            }
          ],
          anyMap: {
            victim: true,
            wound: true
          }
        }
      ]);
      builder.addColumns([stamp('monster', 'name', [])]);

      expect(builder.query()).toEqual([
        'monster',
        [
          'victim',
          ['wound', ['location', '::equals', 'arm'], '::any'],
          '::any'
        ],
        '::all',
        ['name']
      ]);
    });

    it('shallow paths in filters with branching models', () => {
      builder.addRootModel('monster');
      builder.addRecordFilters([
        {
          modelName: 'victim',
          clauses: [
            {
              subclauses: [
                {
                  attributeName: 'name',
                  operator: '::equals',
                  operand: 'Hercules',
                  attributeType: 'text'
                }
              ],
              modelName: 'victim',
              any: true
            }
          ],
          anyMap: {
            victim: true
          }
        }
      ]);
      builder.addColumns([stamp('monster', 'name', [])]);

      expect(builder.query()).toEqual([
        'monster',
        ['victim', ['name', '::equals', 'Hercules'], '::any'],
        '::all',
        ['name']
      ]);
    });

    it('paths that go up and down the tree', () => {
      builder.addRootModel('prize');
      builder.addRecordFilters([
        {
          modelName: 'wound',
          clauses: [
            {
              subclauses: [
                {
                  attributeName: 'location',
                  operator: '::equals',
                  operand: 'arm',
                  attributeType: 'text'
                }
              ],
              modelName: 'wound',
              any: true
            }
          ],
          anyMap: {
            victim: true,
            wound: true
          }
        }
      ]);
      builder.addColumns([stamp('prize', 'name', [])]);

      expect(builder.query()).toEqual([
        'prize',
        [
          'labor',
          [
            'monster',
            [
              'victim',
              ['wound', ['location', '::equals', 'arm'], '::any'],
              '::any'
            ],
            '::any'
          ],
          '::any'
        ],
        '::all',
        ['name']
      ]);
    });

    it('paths that go up and down and terminate in a table', () => {
      builder.addRootModel('monster');
      builder.addRecordFilters([
        {
          modelName: 'prize',
          clauses: [
            {
              subclauses: [
                {
                  attributeName: 'name',
                  operator: '::equals',
                  operand: 'Apples',
                  attributeType: 'text'
                }
              ],
              modelName: 'prize',
              any: true
            }
          ],
          anyMap: {
            prize: true
          }
        }
      ]);
      builder.addColumns([stamp('monster', 'name', [])]);

      expect(builder.query()).toEqual([
        'monster',
        ['labor', ['prize', ['name', '::equals', 'Apples'], '::any'], '::any'],
        '::all',
        ['name']
      ]);
    });

    it('paths that go up the tree', () => {
      builder.addRootModel('prize');
      builder.addRecordFilters([
        {
          modelName: 'labor',
          clauses: [
            {
              subclauses: [
                {
                  attributeName: 'name',
                  operator: '::in',
                  operand: 'Lion,Hydra',
                  attributeType: 'text'
                }
              ],
              modelName: 'labor',
              any: true
            }
          ],
          anyMap: {}
        }
      ]);
      builder.addColumns([stamp('prize', 'name', [])]);

      expect(builder.query()).toEqual([
        'prize',
        ['labor', ['name', '::in', ['Lion', 'Hydra']], '::any'],
        '::all',
        ['name']
      ]);
    });
  });

  describe('handles every for', () => {
    it('deep paths in filters with some non-branching models', () => {
      builder.addRootModel('labor');
      builder.addRecordFilters([
        {
          modelName: 'wound',
          clauses: [
            {
              subclauses: [
                {
                  attributeName: 'location',
                  operator: '::equals',
                  operand: 'arm',
                  attributeType: 'text'
                }
              ],
              modelName: 'wound',
              any: true
            }
          ],
          anyMap: {
            victim: false,
            wound: false
          }
        }
      ]);
      builder.addColumns([stamp('labor', 'name', [])]);

      expect(builder.query()).toEqual([
        'labor',
        [
          'monster',
          [
            'victim',
            ['wound', ['location', '::equals', 'arm'], '::every'],
            '::every'
          ],
          '::any'
        ],
        '::all',
        ['name']
      ]);
    });

    it('deep paths in filters with branching models only', () => {
      builder.addRootModel('monster');
      builder.addRecordFilters([
        {
          modelName: 'wound',
          clauses: [
            {
              subclauses: [
                {
                  attributeName: 'location',
                  operator: '::equals',
                  operand: 'arm',
                  attributeType: 'text'
                }
              ],
              modelName: 'wound',
              any: true
            }
          ],
          anyMap: {
            victim: false,
            wound: false
          }
        }
      ]);
      builder.addColumns([stamp('monster', 'name', [])]);

      expect(builder.query()).toEqual([
        'monster',
        [
          'victim',
          ['wound', ['location', '::equals', 'arm'], '::every'],
          '::every'
        ],
        '::all',
        ['name']
      ]);
    });

    it('shallow paths in filters with branching models', () => {
      builder.addRootModel('monster');
      builder.addRecordFilters([
        {
          modelName: 'victim',
          clauses: [
            {
              subclauses: [
                {
                  attributeName: 'name',
                  operator: '::equals',
                  operand: 'Hercules',
                  attributeType: 'text'
                }
              ],
              modelName: 'victim',
              any: true
            }
          ],
          anyMap: {
            victim: false
          }
        }
      ]);
      builder.addColumns([stamp('monster', 'name', [])]);

      expect(builder.query()).toEqual([
        'monster',
        ['victim', ['name', '::equals', 'Hercules'], '::every'],
        '::all',
        ['name']
      ]);
    });

    it('paths that go up and down the tree', () => {
      builder.addRootModel('prize');
      builder.addRecordFilters([
        {
          modelName: 'wound',
          clauses: [
            {
              subclauses: [
                {
                  attributeName: 'location',
                  operator: '::equals',
                  operand: 'arm',
                  attributeType: 'text'
                }
              ],
              modelName: 'wound',
              any: true
            }
          ],
          anyMap: {
            victim: false,
            wound: false
          }
        }
      ]);
      builder.addColumns([stamp('prize', 'name', [])]);

      expect(builder.query()).toEqual([
        'prize',
        [
          'labor',
          [
            'monster',
            [
              'victim',
              ['wound', ['location', '::equals', 'arm'], '::every'],
              '::every'
            ],
            '::any'
          ],
          '::any'
        ],
        '::all',
        ['name']
      ]);
    });

    it('paths that go up and down and terminate in a table', () => {
      builder.addRootModel('monster');
      builder.addRecordFilters([
        {
          modelName: 'prize',
          clauses: [
            {
              subclauses: [
                {
                  attributeName: 'name',
                  operator: '::equals',
                  operand: 'Apples',
                  attributeType: 'text'
                }
              ],
              modelName: 'prize',
              any: true
            }
          ],
          anyMap: {
            prize: false
          }
        }
      ]);
      builder.addColumns([stamp('monster', 'name', [])]);

      expect(builder.query()).toEqual([
        'monster',
        [
          'labor',
          ['prize', ['name', '::equals', 'Apples'], '::every'],
          '::any'
        ],
        '::all',
        ['name']
      ]);
    });

    it('paths that go up the tree', () => {
      builder.addRootModel('prize');
      builder.addRecordFilters([
        {
          modelName: 'labor',
          clauses: [
            {
              subclauses: [
                {
                  attributeName: 'name',
                  operator: '::in',
                  operand: 'Lion,Hydra',
                  attributeType: 'text'
                }
              ],
              modelName: 'labor',
              any: true
            }
          ],
          anyMap: {}
        }
      ]);
      builder.addColumns([stamp('prize', 'name', [])]);

      expect(builder.query()).toEqual([
        'prize',
        ['labor', ['name', '::in', ['Lion', 'Hydra']], '::any'],
        '::all',
        ['name']
      ]);
    });

    it('correctly handles numeric filters', () => {
      builder.addRootModel('monster');
      builder.addColumns([stamp('monster', 'name', [])]);
      builder.addRecordFilters([
        {
          modelName: 'labor',
          clauses: [
            {
              subclauses: [
                {
                  attributeName: 'number',
                  operator: '::equals',
                  operand: '2',
                  attributeType: 'number'
                }
              ],
              modelName: 'labor',
              any: true
            }
          ],
          anyMap: {}
        },
        {
          modelName: 'labor',
          clauses: [
            {
              subclauses: [
                {
                  attributeName: 'number',
                  operator: '::in',
                  operand: '1,3,5',
                  attributeType: 'number'
                }
              ],
              modelName: 'labor',
              any: true
            }
          ],
          anyMap: {}
        },
        {
          modelName: 'labor',
          clauses: [
            {
              subclauses: [
                {
                  attributeName: 'number',
                  operator: '::notin',
                  operand: '2,4,6',
                  attributeType: 'number'
                }
              ],
              modelName: 'labor',
              any: true
            }
          ],
          anyMap: {}
        },
        {
          modelName: 'labor',
          clauses: [
            {
              subclauses: [
                {
                  attributeName: 'number',
                  operator: '::>=',
                  operand: '5,6',
                  attributeType: 'number'
                }
              ],
              modelName: 'labor',
              any: true
            }
          ],
          anyMap: {}
        }
      ]);

      expect(builder.query()).toEqual([
        'monster',
        [
          '::and',
          ['labor', ['number', '::equals', 2], '::any'],
          ['labor', ['number', '::in', [1, 3, 5]], '::any'],
          ['labor', ['number', '::notin', [2, 4, 6]], '::any'],
          ['labor', ['number', '::>=', 5], '::any']
        ],
        '::all',
        ['name']
      ]);
    });

    it('correctly handles non-numeric filters', () => {
      builder.addRootModel('monster');
      builder.addColumns([stamp('monster', 'name', [])]);
      builder.addRecordFilters([
        {
          modelName: 'labor',
          clauses: [
            {
              subclauses: [
                {
                  attributeName: 'name',
                  operator: '::equals',
                  operand: '2',
                  attributeType: 'text'
                }
              ],
              modelName: 'labor',
              any: true
            }
          ],
          anyMap: {}
        },
        {
          modelName: 'labor',
          clauses: [
            {
              subclauses: [
                {
                  attributeName: 'name',
                  operator: '::in',
                  operand: '1,3,5',
                  attributeType: 'text'
                }
              ],
              modelName: 'labor',
              any: true
            }
          ],
          anyMap: {}
        },
        {
          modelName: 'labor',
          clauses: [
            {
              subclauses: [
                {
                  attributeName: 'name',
                  operator: '::notin',
                  operand: '2,4,6',
                  attributeType: 'text'
                }
              ],
              modelName: 'labor',
              any: true
            }
          ],
          anyMap: {}
        },
        {
          modelName: 'labor',
          clauses: [
            {
              subclauses: [
                {
                  attributeName: 'name',
                  operator: '::>=',
                  operand: '5,6',
                  attributeType: 'text'
                }
              ],
              modelName: 'labor',
              any: true
            }
          ],
          anyMap: {}
        }
      ]);

      expect(builder.query()).toEqual([
        'monster',
        [
          '::and',
          ['labor', ['name', '::equals', '2'], '::any'],
          ['labor', ['name', '::in', ['1', '3', '5']], '::any'],
          ['labor', ['name', '::notin', ['2', '4', '6']], '::any'],
          ['labor', ['name', '::>=', '5,6'], '::any']
        ],
        '::all',
        ['name']
      ]);
    });
  });
});
