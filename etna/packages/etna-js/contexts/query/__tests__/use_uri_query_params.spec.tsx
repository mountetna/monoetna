import React from 'react';

import {render, screen, waitFor, fireEvent} from '@testing-library/react';
import '@testing-library/jest-dom/extend-expect';

import TestComponent from './mock_component';

describe('useUriQueryParams', () => {
  let setColumnsMock: any;
  let setWhereMock: any;
  let setRootMock: any;

  beforeEach(() => {
    setColumnsMock = jest.fn();
    setWhereMock = jest.fn();
    setRootMock = jest.fn();
  });

  describe('old clause format', () => {
    it('correctly migrates', async () => {
      const input = {
        columns: [
          {
            anyMap: {},
            modelName: 'labor',
            slices: [
              {
                clause: {
                  modelName: 'labor',
                  any: true,
                  attributeName: 'name',
                  attributeType: 'string',
                  operator: '::equals',
                  operand: 'foo'
                }
              }
            ]
          }
        ],
        recordFilters: [
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
              },
              {
                modelName: 'labor',
                any: true,
                attributeName: 'name',
                attributeType: 'string',
                operator: '::matches',
                operand: 'foo'
              }
            ]
          }
        ],
        orRecordFilterIndices: [1]
      };

      const {asFragment} = render(
        <TestComponent
          setWhereState={setWhereMock}
          setQueryColumns={setColumnsMock}
          setRootModel={setRootMock}
          input={input}
        />
      );

      await waitFor(() => screen.getByText('Rendered'));

      expect(setWhereMock).toHaveBeenCalledWith({
        recordFilters: [
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
              },
              {
                modelName: 'labor',
                any: true,
                subclauses: [
                  {
                    attributeName: 'name',
                    attributeType: 'string',
                    operator: '::matches',
                    operand: 'foo'
                  }
                ]
              }
            ]
          }
        ],
        orRecordFilterIndices: [1]
      });
      expect(setColumnsMock).toHaveBeenCalledWith([
        {
          anyMap: {},
          modelName: 'labor',
          slices: [
            {
              clause: {
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
            }
          ]
        }
      ]);
    });
  });

  it('generates empty subclause with empty old format', async () => {
    const input = {
      columns: [
        {
          anyMap: {},
          modelName: 'labor',
          slices: [
            {
              clause: {
                modelName: 'labor',
                any: true,
                attributeName: '',
                attributeType: '',
                operator: '',
                operand: ''
              }
            }
          ]
        }
      ],
      recordFilters: [
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
      ],
      orRecordFilterIndices: [1]
    };

    const {asFragment} = render(
      <TestComponent
        setWhereState={setWhereMock}
        setQueryColumns={setColumnsMock}
        setRootModel={setRootMock}
        input={input}
      />
    );

    await waitFor(() => screen.getByText('Rendered'));

    expect(setWhereMock).toHaveBeenCalledWith({
      recordFilters: [
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
      ],
      orRecordFilterIndices: [1]
    });
    expect(setColumnsMock).toHaveBeenCalledWith([
      {
        anyMap: {},
        modelName: 'labor',
        slices: [
          {
            clause: {
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
          }
        ]
      }
    ]);
  });

  it('ignores new clause format', async () => {
    const input = {
      columns: [
        {
          anyMap: {},
          modelName: 'labor',
          slices: [
            {
              clause: {
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
            }
          ]
        }
      ],
      recordFilters: [
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
                },
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
      ],
      orRecordFilterIndices: [1]
    };

    const {asFragment} = render(
      <TestComponent
        setWhereState={setWhereMock}
        setQueryColumns={setColumnsMock}
        setRootModel={setRootMock}
        input={input}
      />
    );

    await waitFor(() => screen.getByText('Rendered'));

    expect(setWhereMock).toHaveBeenCalledWith({
      recordFilters: [
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
                },
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
      ],
      orRecordFilterIndices: [1]
    });
    expect(setColumnsMock).toHaveBeenCalledWith([
      {
        anyMap: {},
        modelName: 'labor',
        slices: [
          {
            clause: {
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
          }
        ]
      }
    ]);
  });
});
