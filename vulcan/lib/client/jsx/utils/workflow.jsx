import React from 'react';

import * as _ from 'lodash';

import XYPlotModel from '../models/xy_plot';

import ListInput from 'etna-js/components/inputs/list_input';
import DropdownInput from 'etna-js/components/inputs/dropdown_input';
import {
  IntegerInput,
  FloatInput
} from 'etna-js/components/inputs/numeric_input';
import SlowTextInput from 'etna-js/components/inputs/slow_text_input';

import {TYPE} from '../models/steps';

export const wrapPaneItem = (item) => {
  return (
    <div className='view_item'>
      <div className='item_name'>{item.name}</div>
      <div className='item_view'>{item.value}</div>
    </div>
  );
};

export const wrapEditableInputs = (inputs, handleInputChange) => {
  return Object.keys(inputs).map((inputName, index) => {
    let input = inputs[inputName];
    let name = input.label || inputName;
    let key = `${inputName}-${index}`;

    switch (input.type) {
      case TYPE.INTEGER:
        return wrapPaneItem({
          name,
          value: (
            <IntegerInput
              key={key}
              defaultValue={input.default}
              onChange={(e) => {
                handleInputChange(inputName, e);
              }}
            ></IntegerInput>
          )
        });
      case TYPE.FLOAT:
        return wrapPaneItem({
          name,
          value: (
            <FloatInput
              key={key}
              defaultValue={input.default}
              onChange={(e) => {
                handleInputChange(inputName, e);
              }}
            ></FloatInput>
          )
        });
      case TYPE.BOOL:
        return wrapPaneItem({
          name,
          value: (
            <input
              key={key}
              type='checkbox'
              className='text_box'
              onChange={(e) => {
                handleInputChange(inputName, e);
              }}
              defaultChecked={input.default}
            />
          )
        });
      case TYPE.MULTISELECT_STRING:
        return wrapPaneItem({
          name,
          value: (
            <ListInput
              key={key}
              placeholder='Select items from the list'
              className='link_text'
              values={input.default || []}
              itemInput={DropdownInput}
              list={input.options || []}
              onChange={(e) => {
                handleInputChange(inputName, e);
              }}
            />
          )
        });
      default:
        return wrapPaneItem({
          name,
          value: (
            <SlowTextInput
              key={key}
              defaultValue={input.default}
              onChange={(e) => {
                handleInputChange(inputName, e);
              }}
            ></SlowTextInput>
          )
        });
    }
  });
};

export const uiStepInputNames = (step) => {
  return step.out.map((output) => `${step.name}/${output}`);
};

export const missingUiInputs = (step, session) => {
  return uiStepInputNames(step).filter(
    (outputName) => !Object.keys(session.inputs).includes(outputName)
  );
};

export const inputNamesToHashStub = (inputNames) => {
  // Convert a list of input strings to
  //   Hash, where all values are `null`.
  return inputNames.reduce((result, input) => {
    if (!result.hasOwnProperty(input)) {
      result[input] = null;
    }

    return result;
  }, {});
};

export const uiStepType = (step) => {
  return step.run.split('/')[1].replace('.cwl', '');
};

export const uiStepOptions = ({step, pathIndex, status}) => {
  // Pull out any previous step's output data that is a required
  //   input into this UI step.
  // Assume a single input data file for now.
  let previousStepName = step.in[0].source[0];
  let outputKey = step.in[0].source[1];

  let previousStep = status[pathIndex].find((s) => s.name === previousStepName);

  if (!previousStep || !previousStep.data || !previousStep.data[outputKey])
    return [];

  let outputData = previousStep.data[outputKey];

  return Array.isArray(outputData) ? outputData : [outputData];
};

export const validStep = ({workflow, pathIndex, stepIndex}) => {
  return !!(
    workflow.steps &&
    workflow.steps[pathIndex] &&
    workflow.steps[pathIndex][stepIndex]
  );
};

export const validPath = ({workflow, pathIndex}) => {
  return !!(workflow.steps && workflow.steps[pathIndex]);
};

export const inputType = ({workflow, input}) => {
  // First check the input itself for a locally typed input.
  // If no type present, check the workflow inputs for
  //   for type information.
};

export const hasUiInput = (step) => {
  return step.run && step.run.startsWith('ui-queries/');
};

export const defaultInputValues = (workflow) => {
  return Object.keys(workflow.inputs).reduce((result, inputName) => {
    if (workflow.inputs[inputName].default) {
      result[inputName] = workflow.inputs[inputName].default;
    }
    return result;
  }, {});
};

export const allInputsDefined = (workflow, userInputs) => {
  // Session Inputs can also include UI interaction inputs,
  //   so they won't be defined in the workflow inputs.
  return (
    Object.keys(workflow.inputs).every((primaryInput) =>
      Object.keys(userInputs).includes(primaryInput)
    ) &&
    Object.keys(userInputs).every((key) => {
      let userInput = userInputs[key];
      // The userInput !== userInput is to check for NaN, because
      //    NaN !== NaN
      // For multiselect, need to make sure the inputs are not
      //    [""] or []...
      return !(
        null === userInput ||
        undefined === userInput ||
        userInput !== userInput ||
        _.isEqual([''], userInput) ||
        _.isEqual([], userInput)
      );
    })
  );
};

export const stepDataUrls = ({workflow, pathIndex, stepIndex}) => {
  return workflow.steps[pathIndex][stepIndex].out.map((s) => s.data_url);
};

export const stepOutputs = (workflow, pathIndex, stepIndex) => {
  if (
    !workflow.steps ||
    !workflow.steps[pathIndex] ||
    !workflow.steps[pathIndex][stepIndex]
  )
    return [];

  return workflow.steps[pathIndex][stepIndex].out;
};

const plotType = (step) => {
  return step.run.replace(/\.cwl/g, '');
};

export const plotModelForStep = (step, consignment, parentWidth) => {
  // Take a consignment and step, and generate
  //   the plot config and data objects
  //   needed to generate a plot.
  // This should mostly be reverse-engineering the
  //   series, data, and label values from the step +
  //   consignment.
  const type = plotType(step);

  let plot;
  if ('xy' === type) {
    plot = new XYPlotModel(step, consignment, parentWidth);
  }

  return plot;
};
