import {VulcanState} from "../reducers/vulcan_reducer";
import {Dispatch, useEffect, useMemo} from "react";
import {patchInputs, removeInputs, VulcanAction} from "../actions/vulcan";
import {
  findSourceDependencies,
  inputValueNonEmpty, isPendingUiQuery,
  missingUiQueryOutputs,
  pendingSteps, sourceNameOfReference, statusOfStep, stepInputDataRaw, stepOfSource, stepOfStatus, uiQueryOfStep
} from "../selectors/workflow_selectors";

export function useInputStateManagement(state: VulcanState, dispatch: Dispatch<VulcanAction>) {
  const {inputs, workflow, status, data, session} = state;

  useEffect(() => {
    if (workflow == null) return;
    const deletes: {[k: string]: true} = {};

    const droppedSteps = workflow.steps[0].filter(step => {
      const stepStatus = statusOfStep(step, status);
      if (step.out.length === 0) return false;

      if (!stepStatus) return true;

      return !step.out.every(outName => {
        const source = sourceNameOfReference([step.name, outName]);
        if (source in inputs) {
          return true;
        }

        if (uiQueryOfStep(step)) return true;

        console.log({step});
        if (!stepStatus.downloads) return false;

        return outName in stepStatus.downloads;
      });
    })

    droppedSteps.forEach(step => {
      step.out.forEach(outName => {
        const source = sourceNameOfReference([step.name, outName]);
        if (source in inputs) deletes[source] = true;
        workflow.dependencies_of_outputs[source].forEach(dependent => {
          if (dependent in inputs) deletes[dependent] = true;
        });
      })
    })

    const d = Object.keys(deletes);
    if (d.length > 0) dispatch(removeInputs(d));
  }, [inputs, workflow, status, data, session]);

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