import {migrateSubclauses} from '../query_uri_params';

describe('migrateSubclauses', () => {
  describe('old clause format', () => {
    it('correctly migrates', () => {
      const input = [
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

      expect(migrateSubclauses(input)).toEqual([
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

    it('correctly migrates with no operand', () => {
      const input = [
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

      expect(migrateSubclauses(input)).toEqual([
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

  it('generates empty subclause with empty old format', () => {
    const input = [
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

    expect(migrateSubclauses(input)).toEqual([
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

  it('ignores new clause format', () => {
    const input = [
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

    expect(migrateSubclauses(input)).toEqual([
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
