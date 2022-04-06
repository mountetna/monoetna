import React, {ReactNode, ReactPortal} from 'react';
import {
  integrateElement,
  setupBefore
} from '../../../../../test_utils/integration';

import {Maybe, some} from '../../../../../selectors/maybe';
import DataTransformationInput, {
  dataFrameJsonToNestedArray,
  nestedArrayToDataFrameJson
} from '../data_transformation';
import {DataEnvelope} from '../input_types';
import {
  clickNode,
  findAllByClassName,
  matchesTextPredicate,
  matchesTypePredicate,
  text
} from '../../../../../test_utils/rendered';
import {ReactTestInstance} from 'react-test-renderer';
import ReactDOM from 'react-dom';

describe('dataFrameJsonToNestedArray', () => {
  it('can reformat the data', () => {
    const inputDF = {
      col1: {
        '0': 1,
        '1': 2,
        '2': 10
      },
      col2: {
        '0': 'abc',
        '1': '21',
        '2': 'xyz'
      }
    };
    const output = dataFrameJsonToNestedArray(some(inputDF));

    expect(output).toEqual([
      ['col1', 'col2'],
      [1, 'abc'],
      [2, '21'],
      [10, 'xyz']
    ]);
  });
});

describe('nestedArrayToDataFrameJson', () => {
  it('can reformat the data', () => {
    const input = [
      ['col1', 'col2'],
      [1, 'abc'],
      [2, '21'],
      [10, 'xyz']
    ];

    const output = nestedArrayToDataFrameJson(input);

    expect(output).toEqual({
      col1: {
        '0': 1,
        '1': 2,
        '2': 10
      },
      col2: {
        '0': 'abc',
        '1': '21',
        '2': 'xyz'
      }
    });
  });
});

describe('DataTransformationInput', () => {
  const onChange = setupBefore(() => jest.fn());
  const value = setupBefore(
    () => null as Maybe<DataEnvelope<{[key: string]: any}>>
  );
  const data = setupBefore(() => {
    return {
      data_frame: {
        record_name_01: {
          '0': 1,
          '1': 0.25,
          '2': 'data'
        },
        record_name_02: {
          '0': 200,
          '1': 1.111,
          '2': '=IF(A2>100, 1, 0)'
        }
      }
    } as DataEnvelope<DataEnvelope<{[key: string]: any}>>;
  });

  const integrated = setupBefore(() =>
    integrateElement(
      <DataTransformationInput
        onChange={onChange.value}
        value={value.value}
        data={data.value}
        numOutputs={2}
      />
    )
  );

  async function clickButton(component: ReactTestInstance, idx: number) {
    await clickNode(component, matchesTypePredicate('button'), idx);
  }

  it('includes a brief description of the data frame size', async () => {
    const {node} = integrated.value;
    expect(
      matchesTextPredicate(
        'Your data frame has 3 rows and 2 columns. You can preview or edit the data frame now, or just click "Commit" to accept the raw data.Review or edit data frame'
      )(node.root)
    ).toEqual(true);
  });

  it('correctly initializes the default value as nested calculated_data / formulaic_data hash', async () => {
    const {node} = integrated.value;
    expect(onChange.value).toHaveBeenCalledWith(
      some({
        calculated_data: some({
          record_name_01: {
            '0': 1,
            '1': 0.25,
            '2': 'data'
          },
          record_name_02: {
            '0': 200,
            '1': 1.111,
            '2': '=IF(A2>100, 1, 0)'
          }
        }),
        formulaic_data: some({
          record_name_01: {
            '0': 1,
            '1': 0.25,
            '2': 'data'
          },
          record_name_02: {
            '0': 200,
            '1': 1.111,
            '2': '=IF(A2>100, 1, 0)'
          }
        })
      }),
      true
    );
  });

  describe('with a buffered input.value', () => {
    value.replace(() => {
      return some({
        calculated_data: some({
          record_name_01: {
            '0': 1,
            '1': 0.25,
            '2': 'data',
            '3': 'abc1232'
          },
          record_name_02: {
            '0': 200,
            '1': 1.111,
            '2': '=IF(A2>100, 1, 0)',
            '3': '-9.9'
          }
        }),
        formulaic_data: some({
          record_name_01: {
            '0': 1,
            '1': 0.25,
            '2': 'data',
            '3': 2
          },
          record_name_02: {
            '0': 200,
            '1': 1.111,
            '2': '=IF(A2>100, 1, 0)',
            '3': '=1+1'
          }
        })
      }) as Maybe<DataEnvelope<{[key: string]: any}>>;
    });

    it('correctly renders', async () => {
      const {node} = integrated.value;
      expect(
        matchesTextPredicate(
          'Your data frame has 4 rows and 2 columns. You can preview or edit the data frame now, or just click "Commit" to accept the raw data.Review or edit data frame** You have modified the data frame. **Revert to raw data'
        )(node.root)
      ).toEqual(true);
    });

    describe('with differently-named CWL outputs', () => {
      value.replace(() => {
        return some({
          a: some({
            record_name_01: {
              '0': 1,
              '1': 0.25,
              '2': 'data',
              '3': 'abc1232',
              '4': '000'
            },
            record_name_02: {
              '0': 200,
              '1': 1.111,
              '2': '=IF(A2>100, 1, 0)',
              '3': '-9.9',
              '4': '111'
            }
          }),
          b: some({
            record_name_01: {
              '0': 1,
              '1': 0.25,
              '2': 'data',
              '3': 2,
              '4': 'abc'
            },
            record_name_02: {
              '0': 200,
              '1': 1.111,
              '2': '=IF(A2>100, 1, 0)',
              '3': '=1+1',
              '4': '=1*9'
            }
          })
        }) as Maybe<DataEnvelope<{[key: string]: any}>>;
      });

      it('correctly renders', async () => {
        const {node} = integrated.value;

        expect(
          matchesTextPredicate(
            'Your data frame has 5 rows and 2 columns. You can preview or edit the data frame now, or just click "Commit" to accept the raw data.Review or edit data frame** You have modified the data frame. **Revert to raw data'
          )(node.root)
        ).toEqual(true);
      });
    });
  });

  describe('with a committed input.value', () => {
    value.replace(() => {
      return some({
        calculated_data: {
          record_name_01: {
            '0': 1,
            '1': 0.25,
            '2': 'data',
            '3': 'abc1232'
          },
          record_name_02: {
            '0': 200,
            '1': 1.111,
            '2': '=IF(A2>100, 1, 0)',
            '3': '-9.9'
          }
        },
        formulaic_data: {
          record_name_01: {
            '0': 1,
            '1': 0.25,
            '2': 'data',
            '3': 2
          },
          record_name_02: {
            '0': 200,
            '1': 1.111,
            '2': '=IF(A2>100, 1, 0)',
            '3': '=1+1'
          }
        }
      }) as Maybe<DataEnvelope<{[key: string]: any}>>;
    });

    it('correctly renders', async () => {
      const {node} = integrated.value;

      expect(
        matchesTextPredicate(
          'Your data frame has 4 rows and 2 columns. You can preview or edit the data frame now, or just click "Commit" to accept the raw data.Review or edit data frame** You have modified the data frame. **Revert to raw data'
        )(node.root)
      ).toEqual(true);
    });

    describe('with differently-named CWL outputs', () => {
      value.replace(() => {
        return some({
          a: {
            record_name_01: {
              '0': 1,
              '1': 0.25,
              '2': 'data',
              '3': 'abc1232',
              '4': '000'
            },
            record_name_02: {
              '0': 200,
              '1': 1.111,
              '2': '=IF(A2>100, 1, 0)',
              '3': '-9.9',
              '4': '111'
            }
          },
          b: {
            record_name_01: {
              '0': 1,
              '1': 0.25,
              '2': 'data',
              '3': 2,
              '4': 'abc'
            },
            record_name_02: {
              '0': 200,
              '1': 1.111,
              '2': '=IF(A2>100, 1, 0)',
              '3': '=1+1',
              '4': '=1*9'
            }
          }
        }) as Maybe<DataEnvelope<{[key: string]: any}>>;
      });

      it('correctly renders', async () => {
        const {node} = integrated.value;

        expect(
          matchesTextPredicate(
            'Your data frame has 5 rows and 2 columns. You can preview or edit the data frame now, or just click "Commit" to accept the raw data.Review or edit data frame** You have modified the data frame. **Revert to raw data'
          )(node.root)
        ).toEqual(true);
      });
    });
  });

  // Not clear that this test can work, yet, due to
  //   issues with MUI portals and react-test-renderer
  // * https://github.com/facebook/react/issues/11565
  // * https://github.com/mui/material-ui/issues/12237
  xit('opens a dialog with a table in it', async () => {
    const oldCreatePortal = ReactDOM.createPortal;
    ReactDOM.createPortal = (node: ReactNode): ReactPortal =>
      node as ReactPortal;

    const {node} = integrated.value;
    await clickButton(node.root, 0);
    expect(findAllByClassName(node.root, 'handsontable').length > 0).toEqual(
      true
    );

    ReactDOM.createPortal = oldCreatePortal;
  });
});
