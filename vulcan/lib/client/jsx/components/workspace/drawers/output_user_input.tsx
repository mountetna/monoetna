import React, {useContext} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';
import {OUTPUTS} from '../../ui_components';

import StepIconName from './step_elements/step_icon_name';

import {uiComponentOfStep} from '../../../selectors/workflow_selectors';
import {VulcanConfigElement, WorkspaceStep} from '../../../api_types';
import { fillInputData } from '../ui_definitions/input_types';

export default function OutputUI({step}: {step: WorkspaceStep}) {
  const {state} = useContext(VulcanContext);
  const {workspace, status} = state;
  const config = workspace?.vulcan_config ? workspace.vulcan_config : undefined;

  if (!config) return null;

  const uiOutput = uiComponentOfStep(step.name, config);

  if (!uiOutput) return null;
  const stepType = uiOutput in OUTPUTS ? uiOutput : 'default';

  const data = fillInputData(step as VulcanConfigElement, step, status.last_params, status.file_contents);

  let OutputComponent = OUTPUTS[stepType];

  return (
    <div className='step-output'>
      <StepIconName step={step} />
      <div className='outputs-pane'>
        <OutputComponent data={data}/>
      </div>
    </div>
  );
}
