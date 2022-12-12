import {QueryDeserializer} from '../query_deserializer';

describe('QueryDeserializer', () => {
  it('identifies the raw filters', () => {
    const deserializer = new QueryDeserializer(
      [
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
      ],
      [
        'monster.species',
        'monster.stats',
        'labor.year',
        'labor.constributions.Athens',
        'labor.contributions.Sidon',
        'prize.value',
        'victim.country'
      ]
    );

    expect(deserializer.rawFilters()).toEqual([
      '::and',
      ['labor', ['name', '::in', ['lion', 'hydra', 'apples']], '::any'],
      ['name', '::equals', 'Nemean Lion'],
      ['labor', ['number', '::equals', 2], '::any'],
      ['labor', ['prize', ['name', '::equals', 'Apples'], '::any'], '::any']
    ]);
  });

  it('identifies the raw columns', () => {
    const deserializer = new QueryDeserializer(
      [
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
      ],
      [
        'monster.species',
        'monster.stats',
        'labor.year',
        'labor.constributions.Athens',
        'labor.contributions.Sidon',
        'prize.value',
        'victim.country'
      ]
    );

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

  describe('orRecordFilterIndices', () => {
    it('gives them correct indices when no ::and filters', () => {
      const deserializer = new QueryDeserializer(
        [
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
        ],
        [
          'monster.species',
          'monster.stats',
          'labor.year',
          'labor.constributions.Athens',
          'labor.contributions.Sidon',
          'prize.value',
          'victim.country'
        ]
      );

      expect(deserializer.orRecordFilterIndices()).toEqual([0, 1]);
    });

    it('gives them correct indices when at the end of filter list', () => {
      const deserializer = new QueryDeserializer(
        [
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
        ],
        [
          'monster.species',
          'monster.stats',
          'labor.year',
          'labor.constributions.Athens',
          'labor.contributions.Sidon',
          'prize.value',
          'victim.country'
        ]
      );

      expect(deserializer.orRecordFilterIndices()).toEqual([2, 3]);
    });
  });

  it('identifies rootModel', () => {
    const deserializer = new QueryDeserializer(
      [
        'monster',
        [
          '::and',
          ['labor', ['name', '::in', ['lion', 'hydra', 'apples']], '::any'],
          ['name', '::equals', 'Nemean Lion']
        ],
        '::all',
        ['name']
      ],
      ['monster.name']
    );

    expect(deserializer.rootModel()).toEqual('monster');
  });

  describe('#flattenFilter', () => {
    it('handles nested ::and and ::or filters', () => {
      const deserializer = new QueryDeserializer(
        [
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
        ],
        ['monster.name']
      );

      expect(
        deserializer.flattenFilters([
          '::and',
          ['labor', ['name', '::in', ['lion', 'hydra', 'apples']], '::any'],
          ['name', '::equals', 'Nemean Lion'],
          [
            '::or',
            ['labor', ['year', '::<', 2022], '::every'],
            ['labor', ['number', '::equals', 2], '::any']
          ]
        ])
      ).toEqual([
        {
          anyMap: {
            labor: true
          },
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
              any: false,
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
            labor: true
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

  describe('recordFilters', () => {
    it('works with just one on root model', () => {
      const deserializer = new QueryDeserializer(
        ['monster', ['name', '::equals', 'Nemean Lion'], '::all', ['name']],
        ['monster.name']
      );

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

    it('works with multiple model filters joined by ::and', () => {});

    it('correctly sets ::every in anyMap', () => {});

    it('includes multiple subclauses in root model filter', () => {});

    it('includes multiple subclauses in non-root model filter', () => {});

    it('comma-joins ::in operands', () => {
      const deserializer = new QueryDeserializer(
        [
          'monster',
          ['name', '::in', ['Nemean Lion', 'Lernean Hydra']],
          '::all',
          ['name']
        ],
        ['monster.name']
      );

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

    it('correctly constructs inverted operators', () => {});

    it('correctly constructs terminal operators', () => {});
  });

  describe('columns', () => {
    it('works with no edited user_columns', () => {});

    it('works with edited user_columns', () => {});

    it('comma-joins ::slice operands', () => {});
  });
});
