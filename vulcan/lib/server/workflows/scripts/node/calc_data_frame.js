import * as dataflow from 'archimedes-node/functions/dataflow/index.js';
import * as dfUtils from 'archimedes-node/functions/dataframe/index.js';

// const originalData = await dataflow.inputJson('original_data.json');
// const userData = await dataflow.inputJson('user_data.json');

let mockOrig = {
  col1: {
    '0': 1,
    '1': 2,
    '2': 3
  },
  col2: {
    '0': 'a',
    '1': 'b',
    '2': 'c'
  }
};

let mockUser = {
  col1: {
    '0': 1,
    '1': 2,
    '2': 3
  },
  col2: {
    '0': 'a',
    '1': 'b',
    '2': 'c'
  },
  col3: {
    '0': '=A2+5',
    '1': '=A3+5',
    '2': ''
  }
};

dataflow.outputJson(
  dfUtils.zipDF(mockOrig, mockUser),
  'full_user_data.json',
  {},
  '/app/'
);

// const originalData = await dataflow.inputJson('original_data.json');
// const userData = await dataflow.inputJson('user_data.json');

// dataflow.outputJson(
//   dfUtils.zipDF(originalData, userData),
//   'full_user_data.json'
// );
