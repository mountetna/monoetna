import React from 'react';
import {StylesProvider, StylesOptions} from '@material-ui/styles/';
import {FlagsProvider, FlagsState} from '../flags/flags-context';

const generateClassName: StylesOptions['generateClassName'] = (
  rule,
  sheet
): string => `${sheet!.options.classNamePrefix}-${rule.key}`;

export const flagsSpecWrapper =
  (mockState: FlagsState) =>
  ({children}: {children?: any}) =>
    (
      <StylesProvider generateClassName={generateClassName}>
        <FlagsProvider state={mockState}>{children}</FlagsProvider>
      </StylesProvider>
    );
