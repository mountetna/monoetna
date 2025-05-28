import React, {ReactNode, ReactPortal} from 'react';
import {some} from '../../../../selectors/maybe';
import {fillInputData} from '../input_types';
import { test_step_vc } from '../../../../test_utils/fixtures/workspaces-response';
import { WorkspaceStep } from '../../../../api_types';

describe('fillInputData', () => {
  describe('fills data with ', () => {
    const conf = test_step_vc('dropdown', {files: ['a.txt']}, {files: ['out.txt']})
    const params = {a: ['a', 'b', 'c']}
    const files = {'a.txt': ['a', 'b', 'c']}
    fillInputData(
      conf, conf as WorkspaceStep,
      params,
      files
    )
  });
});
