#!/usr/bin/env node

import * as dataflow from 'archimedes-node/functions/dataflow/index.js';
import * as dfUtils from 'archimedes-node/functions/dataframe/index.js';

const originalData = await dataflow.inputJson('original_data.json');
const userData = await dataflow.inputJson('user_data.json');

dataflow.outputJson(
  dfUtils.zipDF(originalData, userData),
  'full_user_data.json'
);
