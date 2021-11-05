import React from 'react';
import thunk from 'redux-thunk';
import configureMockStore from 'redux-mock-store';
import {Provider} from 'react-redux';
import {StylesOptions, StylesProvider} from '@material-ui/styles/';
import {Store} from 'redux';

import {Model} from 'etna-js/models/magma-model';
import {MagmaProvider} from 'etna-js/contexts/magma-context';

export const mockStore = configureMockStore([thunk]);

export const mockDate = () => {
  const currentDate = new Date();

  global.Date = jest.fn(() => currentDate) as any;
};

export const mockFetch = () => (global.fetch = fetch);

export const generateClassName: StylesOptions['generateClassName'] = (
  rule,
  sheet
): string => `${sheet!.options.classNamePrefix}-${rule.key}`;

export const redcapSpecWrapper =
  ({props}: {props: {[key: string]: any}}) =>
  ({children}: {children?: any}) =>
    (
      <StylesProvider generateClassName={generateClassName}>
        <MagmaProvider {...props}>{children}</MagmaProvider>
      </StylesProvider>
    );
