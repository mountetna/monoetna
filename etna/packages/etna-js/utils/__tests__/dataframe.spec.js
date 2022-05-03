import {zipJsonDF} from '../dataframe';

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
