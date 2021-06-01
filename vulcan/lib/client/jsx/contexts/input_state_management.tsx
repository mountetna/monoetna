import {VulcanState} from '../reducers/vulcan_reducer';
import {Dispatch, useEffect} from 'react';
import {removeDownloads, removeInputs, VulcanAction} from '../actions/vulcan';
import {
  sourceNameOfReference,
  statusOfStep,
  stepOfSource
} from '../selectors/workflow_selectors';

export function useInputStateManagement(
  state: VulcanState,
  dispatch: Dispatch<VulcanAction>,
  statusIsFresh: boolean
) {
  const {inputs, workflow, status, data, session} = state;
  const {inputs: sessionInputs} = session;

  useEffect(() => {
    console.log('in effect', statusIsFresh, inputs);
    if (workflow == null) return;
    if (!statusIsFresh) return;
    const inputDeletes: {[k: string]: true} = {};
    const downloadDeletes: {[k: string]: true} = {};
    const stepsWithDownloads: {[k: string]: true} = {};

    const droppedSteps = workflow.steps[0].filter((step) => {
      const stepStatus = statusOfStep(step, status);
      if (step.out.length === 0) return false;

      if (!stepStatus) return true;

      return !step.out.every((outName) => {
        const source = sourceNameOfReference([step.name, outName]);
        if (source in sessionInputs) {
          return true;
        }

        if (!stepStatus.downloads) return false;
        stepsWithDownloads[step.name] = true;

        return outName in stepStatus.downloads;
      });
    });

    console.log('droppedSteps', droppedSteps);

    droppedSteps.forEach((step) => {
      step.out.forEach((outName) => {
        const source = sourceNameOfReference([step.name, outName]);
        if (source in sessionInputs) inputDeletes[source] = true;
        if (step.name in stepsWithDownloads) downloadDeletes[step.name] = true;

        workflow.dependencies_of_outputs[source].forEach((dependent) => {
          const stepName = stepOfSource(dependent);
          if (dependent in inputs) inputDeletes[dependent] = true;
          if (stepName && stepName in stepsWithDownloads)
            downloadDeletes[stepName] = true;
        });
      });
    });

    let d = Object.keys(inputDeletes);
    if (d.length > 0) dispatch(removeInputs(d));

    d = Object.keys(downloadDeletes);
    if (d.length > 0) dispatch(removeDownloads(d));
  }, [inputs, workflow, status, data, sessionInputs, statusIsFresh, dispatch]);
}
