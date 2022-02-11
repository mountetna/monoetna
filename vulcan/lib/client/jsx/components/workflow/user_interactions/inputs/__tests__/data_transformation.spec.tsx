import React from 'react';
import {toNestedArray, toJson} from '../data_transformation';

describe('toNestedArray', () => {
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
    const output = toNestedArray(inputDF);

    expect(output).toEqual([
      ['col1', 'col2'],
      [1, 'abc'],
      [2, '21'],
      [10, 'xyz']
    ]);
  });
});

describe('toJson', () => {
  it('can reformat the data', () => {
    const input = [
      ['col1', 'col2'],
      [1, 'abc'],
      [2, '21'],
      [10, 'xyz']
    ];

    const output = toJson(input);

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
