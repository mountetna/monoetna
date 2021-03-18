import * as _ from 'lodash';

import XYPlotModel from '../models/xy_plot';

import {RUN, OUTPUT_COMPONENT} from '../models/steps';

export const stringify = (text) => {
  // For previewing data inputs / outputs, we just want a string
  //   representation. Unsure if input is JSON or a string.

  try {
    text = JSON.stringify(text);
  } catch (e) {}

  return text;
};

export const workflowName = (workflow) =>
  workflow && workflow.name ? workflow.name.replace('.cwl', '') : null;

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

export const uiStepInputDataLink = ({step, pathIndex, status}) => {
  // Pull out any previous step's output data link that is a required
  //   input into this UI step.
  // Assume a single input data link for now.
  let previousStepName = step.in[0].source[0];
  let outputKey = step.in[0].source[1];

  let previousStep = status[pathIndex].find((s) => s.name === previousStepName);

  if (
    !previousStep ||
    !previousStep.downloads ||
    !previousStep.downloads[outputKey]
  )
    return null;

  return previousStep.downloads[outputKey];
};

export const uiStepInputDataRaw = ({step, pathIndex, status}) => {
  // Pull out any previous step's output data that is a required
  //   input into this UI step.
  // Assume a single input data file for now.
  // Provide the raw data object back.
  let previousStepName = step.in[0].source[0];
  let outputKey = step.in[0].source[1];

  let previousStep = status[pathIndex].find((s) => s.name === previousStepName);

  if (!previousStep || !previousStep.data || !previousStep.data[outputKey])
    return null;

  return previousStep.data[outputKey];
};

export const uiStepOptions = ({step, pathIndex, status}) => {
  // Pull out any previous step's output data that is a required
  //   input into this UI step.
  // Assume a single input data file for now.
  const rawData = uiStepInputDataRaw({step, pathIndex, status});
  if (null === rawData) return [];

  return Array.isArray(rawData) ? rawData : [rawData];
};

export const validStep = ({workflow, pathIndex, stepIndex}) => {
  return !!(
    workflow.steps &&
    workflow.steps[pathIndex] &&
    workflow.steps[pathIndex][stepIndex]
  );
};

export const inputGroupName = (input) => {
  let groupName = input.name.split('__')[0];
  if (groupName === input.name) groupName = 'Inputs';

  return groupName;
};

export const validPath = ({workflow, pathIndex}) => {
  return !!(workflow.steps && workflow.steps[pathIndex]);
};

export const hasUiInput = (step) => {
  return step.run && step.run.startsWith(RUN.UI_QUERY);
};

export const hasUiOutput = (step) => {
  return step.run && step.run.startsWith(RUN.UI_OUTPUT);
};

export const defaultInputValues = (workflow) => {
  return Object.keys(workflow.inputs).reduce((result, inputName) => {
    if (null != workflow.inputs[inputName].default) {
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
        null == userInput ||
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

const stepDependsOn = (step, otherStep) => {
  return (
    step.in &&
    step.in.filter((input) => {
      return (
        otherStep.name === input.source[0] &&
        otherStep.out[0] === input.source[1]
      );
    }).length > 0
  );
};

export const shouldDownloadStepData = ({workflow, pathIndex, stepIndex}) => {
  // Only download step data if its output is an input to
  //   a UI INPUT widget or a UI OUTPUT step that is not a LINK
  let step = workflow.steps[pathIndex][stepIndex];
  let dependentSteps = workflow.steps[pathIndex].filter((s) => {
    return (
      (hasUiInput(s) ||
        (hasUiOutput(s) &&
          OUTPUT_COMPONENT.LINK !== s.run.split('/')[1].replace('.cwl', ''))) &&
      stepDependsOn(s, step)
    );
  });

  return dependentSteps.length > 0;
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
