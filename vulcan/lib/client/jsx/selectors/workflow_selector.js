import XYPlotModel from '../models/xy_plot';

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

export const hasUiInput = ({input}) => {
  return input.name.startsWith('ui_');
};

export const inputs = ({session}) => {
  return session.inputs;
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
