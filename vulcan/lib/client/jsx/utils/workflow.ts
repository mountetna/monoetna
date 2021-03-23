import * as _ from 'lodash';

import XYPlotModel from '../models/xy_plot';

import {RUN, OUTPUT_COMPONENT} from '../models/steps';
import {SessionStatusResponse, StepStatus, VulcanSession, Workflow, WorkflowStep} from "../api/types";

export const stringify = (text: any) => {
  // For previewing data inputs / outputs, we just want a string
  //   representation. Unsure if input is JSON or a string.

  try {
    text = JSON.stringify(text);
  } catch (e) {}

  return text;
};

export const uiStepInputNames = (step: WorkflowStep) => {
  return step.out.map((output) => `${step.name}/${output}`);
};

export const missingUiInputs = (step: WorkflowStep, session: VulcanSession) => {
  return uiStepInputNames(step).filter(
    (outputName) => !Object.keys(session.inputs).includes(outputName)
  );
};

export const inputNamesToHashStub = (inputNames: string[]) => {
  // Convert a list of input strings to
  //   Hash, where all values are `null`.
  return inputNames.reduce((result, input) => {
    if (!result.hasOwnProperty(input)) {
      result[input] = null;
    }

    return result;
  }, {} as {[k: string]: any});
};

export const uiStepOptions = ({step, pathIndex, status}) => {
  // Pull out any previous step's output data that is a required
  //   input into this UI step.
  // Assume a single input data file for now.
  const rawData = uiStepInputDataRaw({step, pathIndex, status});
  if (null === rawData) return [];

  return Array.isArray(rawData) ? rawData : [rawData];
};


export const inputGroupName = (input) => {
  let groupName = input.name.split('__')[0];
  if (groupName === input.name) groupName = 'Inputs';

  groupName = groupName.replace(/_/g, ' ');

  return groupName;
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
