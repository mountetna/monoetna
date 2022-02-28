import React, {ReactNode, ReactPortal} from 'react';
import {some} from '../../../../../selectors/maybe';

import {collapseInputValues} from '../input_types';

describe('collapseInputValues', () => {
  describe('for buffer inputs', () => {
    const bufferInputs = {
      aBool: [false],
      anInt: [1],
      aBoolWithoutDefault: null,
      anIntWithoutDefault: null,
      'aStep/a': ['blah'],
      'bStep/a': [123],
      'bStep/b': [321]
    };

    it('wraps multiple, grouped inputs with some', () => {
      const result = collapseInputValues('bStep', 'bStep', bufferInputs, true);

      expect(result).toEqual(
        some({
          a: [123],
          b: [321]
        })
      );
    });

    it('does not wrap singleton inputs with some', () => {
      const result = collapseInputValues('aStep', 'aStep', bufferInputs, true);

      expect(result).toEqual(['blah']);
    });

    it('does not group primary inputs with similar names', () => {
      const result = collapseInputValues(
        undefined,
        'aBool',
        bufferInputs,
        true
      );

      expect(result).toEqual([false]);
    });
  });

  describe('for session inputs', () => {
    const sessionInputs = {
      aBool: false,
      anInt: 1,
      aBoolWithoutDefault: null,
      anIntWithoutDefault: null,
      'aStep/a': 'blah',
      'bStep/a': 123,
      'bStep/b': 321
    };

    it('does not wrap multiple, grouped inputs as some', () => {
      const result = collapseInputValues(
        'bStep',
        'bStep',
        sessionInputs,
        false
      );

      expect(result).toEqual({
        a: 123,
        b: 321
      });
    });

    it('does not wrap singleton inputs with some', () => {
      const result = collapseInputValues(
        'aStep',
        'aStep',
        sessionInputs,
        false
      );

      expect(result).toEqual('blah');
    });

    it('does not group primary inputs with similar names', () => {
      const result = collapseInputValues(
        undefined,
        'aBool',
        sessionInputs,
        false
      );

      expect(result).toEqual(false);
    });
  });
});
