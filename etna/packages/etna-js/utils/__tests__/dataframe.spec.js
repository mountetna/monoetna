import {
  zipJsonDF,
  dataFrameJsonToNestedArray,
  nestedArrayToDataFrameJson,
  merge
} from '../dataframe';

import {some} from '../../selectors/maybe';

describe('dataFrameJsonToNestedArray', () => {
  it('can reformat the data', () => {
    const inputDF = {
      col1: {
        0: 1,
        1: 2,
        2: 10
      },
      col2: {
        0: 'abc',
        1: '21',
        2: 'xyz'
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
        0: 1,
        1: 2,
        2: 10
      },
      col2: {
        0: 'abc',
        1: '21',
        2: 'xyz'
      }
    });
  });
});

describe('merge', () => {
  it('correctly adds original and user data frames', () => {
    const mockOrig = [
      ['col1', 'col2'],
      [1, 'a'],
      [2, 'b'],
      [3, 'c'],
      [4, 'd']
    ];

    const mockUser = [
      ['col1', 'col2', 'col3'],
      [1, 'a', '=A1 * 2'],
      [2, 'b', '=A2 * 2'],
      [3, 'c', null]
    ];

    let result = merge({original: mockOrig, user: mockUser});

    expect(result).toEqual([
      ['col1', 'col2', 'col3'],
      [1, 'a', '=A1 * 2'],
      [2, 'b', '=A2 * 2'],
      [3, 'c', null],
      [4, 'd', null]
    ]);
  });
});

describe('DataFrame Utils', () => {
  it('correctly extends and evaluates', () => {
    let mockOrig = {
      col1: {
        0: 1,
        1: 2,
        2: 3
      },
      col2: {
        0: 'a',
        1: 'b',
        2: 'c'
      }
    };

    let mockUser = {
      col1: {
        0: 1,
        1: 2,
        2: 3
      },
      col2: {
        0: 'a',
        1: 'b',
        2: 'c'
      },
      col3: {
        0: '=A2+5',
        1: '=A3+5',
        2: ''
      }
    };

    let result = zipJsonDF({original: mockOrig, user: mockUser});

    expect(result).toEqual({
      col1: {0: 1, 1: 2, 2: 3},
      col2: {0: 'a', 1: 'b', 2: 'c'},
      col3: {0: 6, 1: 7, 2: 8}
    });
  });
});
