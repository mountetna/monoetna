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

export const isMissingDrawerInputs = (step, session) => {
  return (
    step.in.filter((input) => {
      let inputName = input.source.join('/');
      return (
        !Object.keys(session.inputs).includes(inputName) ||
        !isDefined(session.inputs[inputName])
      );
    }).length > 0
  );
};

export const isMissingStandardInputs = (step, session) => {
  return (
    uiStepInputNames(step).filter((outputName) => {
      console.log(session.inputs[outputName]);
      return (
        !Object.keys(session.inputs).includes(outputName) ||
        !isDefined(session.inputs[outputName])
      );
    }).length > 0
  );
};

export const localStorageKey = (workflow) => `${workflow.name}.session`;

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

export const uiInputDataReady = ({step, pathIndex, status}) => {
  // Check if data has been set for a UI Input step
  return step.in.every((input) => {
    let previousStepName = input.source[0];
    let outputKey = input.source[1];

    let previousStep = status[pathIndex].find(
      (s) => s.name === previousStepName
    );

    if (!previousStep || !previousStep.data || !previousStep.data[outputKey])
      return false;

    return previousStep.data[outputKey];
  });
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

  groupName = groupName.replace(/_/g, ' ');

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

const isDefined = (input) => {
  // The userInput !== userInput is to check for NaN, because
  //    NaN !== NaN
  // For multiselect, need to make sure the inputs are not
  //    [""] or []...
  return !(
    null == input ||
    input !== input ||
    _.isEqual([''], input) ||
    _.isEqual([], input)
  );
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
      return isDefined(userInput);
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

const dependentUiInputNames = ({workflow, inputName}) => {
  let dependentInputNames = [];
  // Start iterating from index 1, since 0 is work path.
  for (let i = 1; i < workflow.steps.length; i++) {
    let path = workflow.steps[i];
    let stepIndex = path.findIndex((s) => s.name === inputName);
    if (stepIndex > -1) {
      for (let j = stepIndex + 1; j < path.length; j++) {
        let futureStep = path[j];
        if (hasUiInput(futureStep)) {
          dependentInputNames = dependentInputNames.concat(
            uiStepInputNames(futureStep)
          );
        }

        // Check other paths for downstream effects
        dependentInputNames = dependentInputNames.concat(
          dependentUiInputNames({workflow, inputName: futureStep.name})
        );
      }
    }
  }

  return dependentInputNames;
};

export const dependentInputsData = ({workflow, inputName, status}) => {
  let dependentNames = dependentUiInputNames({workflow, inputName});

  return dependentNames.reduce((result, dependentName) => {
    let step = workflow.steps[0].find(
      (s) => s.name === dependentName.split('/')[0]
    );
    result.push({
      url: uiStepInputDataLink({step, pathIndex: 0, status}),
      data: null
    });
    return result;
  }, []);
};

export const removeDependentInputs = ({workflow, inputName, userInputs}) => {
  // Search the non-work paths for any UI-query steps that
  //   follow the given input. For any UI-query steps
  //   in those paths, return the userInputs without their
  //   input values.
  let validInputs = {...userInputs};
  // Start iterating from index 1, since 0 is work path.
  dependentUiInputNames({workflow, inputName}).forEach((inputName) => {
    delete validInputs[inputName];
  });
  return validInputs;
};

export const groupUiSteps = (uiSteps) => {
  // Takes an Array of UI steps and checks their names.
  // If they start with our "group" escape sequence
  //   (groupName__rest_of_step_name)
  // this method will group them into a new "step"
  //   where the name is the given `groupName` and
  //   inputs are merged into a single step.
  // NOTE: status, docs, and other things will
  //   be assumed to be the first one found by
  //   the code ... so they are assumed in sync.
  return Object.values(
    uiSteps.reduce((result, step) => {
      let groupName = inputGroupName(step.step);
      if (groupName === 'Inputs') {
        // Not in a group
        result[step.step.name] = step;
      } else {
        // Will just have to default to the first label.
        // No way to resolve any conflict between the step labels.
        if (Object.keys(result).indexOf(groupName) === -1) {
          result[groupName] = {
            step: {
              name: groupName,
              isGroup: true,
              label: groupName,
              run: step.step.run,
              in: []
            },
            index: step.index
          };
        }

        let stepInputs = uiStepInputNames(step.step);
        result[groupName].step.in = result[groupName].step.in.concat(
          stepInputs.map((i) => {
            return {
              source: [i.split('/')[0], i.split('/')[1]],
              doc: step.step.doc,
              label: step.step.label || i
            };
          })
        );
      }

      return result;
    }, {})
  );
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
