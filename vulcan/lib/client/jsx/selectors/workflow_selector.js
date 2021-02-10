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

export const plotDataForStep = (step, consignment) => {
  // Take a consignment and step, and generate
  //   the plot config and data objects
  //   needed to generate a plot.
  // This should mostly be reverse-engineering the
  //   series, data, and label values from the step +
  //   consignment.
};

export const plotLayout = (plot, dimensions) => {
  // Calculate the plot layout given the dimensions
  //   of a parent container.
};
