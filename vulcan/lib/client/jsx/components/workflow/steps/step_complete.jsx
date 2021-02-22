import React, {useState, useContext, useEffect} from 'react';

import {VulcanContext} from '../../../contexts/vulcan';

import ListInput from 'etna-js/components/inputs/list_input';
import SelectInput from 'etna-js/components/inputs/select_input';
import Dropdown from 'etna-js/components/inputs/dropdown';
import {
  IntegerInput,
  FloatInput
} from 'etna-js/components/inputs/numeric_input';
import Toggle from 'etna-js/components/inputs/toggle';

import {TYPE} from '../../../models/steps';
import {inputType, validStep} from '../../../selectors/workflow_selector';

export default function StepComplete() {
  const {workflow, pathIndex, stepIndex, session, status} = useContext(
    VulcanContext
  );

  if (!validStep({workflow, pathIndex, stepIndex})) return null;
  if (!session || !status) return null;

  let step = workflow.steps[pathIndex][stepIndex];

  function wrapPaneItem(item) {
    return (
      <div className='view_item'>
        <div className='item_name'>{item.name}</div>
        <div className='item_view'>{item.value}</div>
      </div>
    );
  }

  // Find the step's inputs + values from the Primary Inputs
  let inputValues = step.in.map((input) => {
    let value;
    let name = input.id;
    if ('primary_inputs' === input.source[0]) {
      let primaryInputName = input.source[1];
      value = session.inputs[primaryInputName];
    } else {
      let outputStepName = input.source[0];
      let outputStepIndex = workflow.steps[pathIndex].findIndex(
        (s) => outputStepName === s.name
      );
      let outputVariableName = input.source[1];
      value = status[(pathIndex, outputStepIndex)].data[outputVariableName];
    }

    return {
      name,
      value
    };
  });

  // For the output values, we'll have to figure out the type
  //   of thing to render. For now, just show the raw data.
  let outputValue;
  let outputName = step.out[0];

  if (
    status[pathIndex][stepIndex].data &&
    status[pathIndex][stepIndex].data[outputName]
  ) {
    outputValue = status[pathIndex][stepIndex].data[outputName];
  }

  return (
    <div>
      <div className='step-inputs inputs-pane'>
        <div class='title'>Step inputs</div>
        <div className='step-inputs-container items'>
          {inputValues.map((input) => {
            return wrapPaneItem(input);
          })}
        </div>
      </div>
      <div className='step-outputs outputs-pane'>
        <div class='title'>Step outputs</div>
        <div className='step-outputs-container items'>
          {wrapPaneItem({
            name: outputName,
            value: outputValue
          })}
        </div>
      </div>
    </div>
  );
}
