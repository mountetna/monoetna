import {VulcanState} from "../reducers/vulcan_reducer";
import {Dispatch, useEffect, useMemo} from "react";
import {patchInputs, removeDownloads, removeInputs, VulcanAction} from "../actions/vulcan";
import {
  isPendingUiQuery,
  missingUiQueryOutputs,
  pendingSteps, sourceNameOfReference, statusOfStep, stepOfSource
} from "../selectors/workflow_selectors";

export function useInputStateManagement(state: VulcanState, dispatch: Dispatch<VulcanAction>, statusIsFresh: boolean) {
  const {inputs, workflow, status, data, session} = state;

  useEffect(() => {
    if (workflow == null) return;
    if (!statusIsFresh) return;
    const inputDeletes: {[k: string]: true} = {};
    const downloadDeletes: {[k: string]: true} = {};
    const stepsWithDownloads: {[k: string]: true} = {};

    const droppedSteps = workflow.steps[0].filter(step => {
      const stepStatus = statusOfStep(step, status);
      if (step.out.length === 0) return false;

      if (!stepStatus) return true;

      return !step.out.every(outName => {
        const source = sourceNameOfReference([step.name, outName]);
        if (source in session.inputs) {
          return true;
        }

        if (!stepStatus.downloads) return false;
        stepsWithDownloads[step.name] = true;

        return outName in stepStatus.downloads;
      });
    });

    droppedSteps.forEach(step => {
      step.out.forEach(outName => {
        const source = sourceNameOfReference([step.name, outName]);
        if (source in session.inputs) inputDeletes[source] = true;
        if (step.name in stepsWithDownloads) downloadDeletes[step.name] = true;

        workflow.dependencies_of_outputs[source].forEach(dependent => {
          const stepName = stepOfSource(dependent);
          if (dependent in inputs) inputDeletes[dependent] = true;
          if (stepName && stepName in stepsWithDownloads) downloadDeletes[stepName] = true;
        });
      })
    })

    let d = Object.keys(inputDeletes);
    if (d.length > 0) dispatch(removeInputs(d));

    d = Object.keys(downloadDeletes);
    if (d.length > 0) dispatch(removeDownloads(d));
  }, [inputs, workflow, status, data, session.inputs, statusIsFresh]);

  // We inject a `null` input into the session,
  //   to indicate that we're waiting for a user input value
  //   to return to the server.
  useEffect(() => {
    if (!workflow) return;

    let newInputs = {};
    let nextUiSteps = pendingSteps(workflow, status).filter(({step}) => isPendingUiQuery(step, status, data, session));

    nextUiSteps.forEach((nextInputStep) => {
      let missingInputs = missingUiQueryOutputs(nextInputStep.step, inputs);

      if (Object.keys(missingInputs).length > 0) {
        // Make sure to copy over the current inputs, otherwise
        //   they'll get wiped out in the reducer.
        newInputs = {
          ...newInputs,
          ...missingInputs,
        };
      }
    });

    if (Object.keys(newInputs).length > 0) dispatch(patchInputs(newInputs));
  }, [inputs, workflow, status, data, session]);
}