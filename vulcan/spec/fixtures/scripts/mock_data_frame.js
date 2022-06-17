#!/usr/bin/env node

import * as dataflow from 'archimedes-node/functions/dataflow/index.js';

dataflow.outputJson(
  {
    col1: {
      '0': 1,
      '1': 2
    },
    col2: {
      '0': 'abc',
      '1': 'xyz'
    }
  },
  'df.json'
);
