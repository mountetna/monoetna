import {VulcanState} from "../reducers/vulcan_reducer";
import {Dispatch, useEffect} from "react";
import {patchInputs, removeInputs, VulcanAction} from "../actions/vulcan";
import {
  findSourceDependencies,
  inputValueNonEmpty, isPendingUiQuery,
  missingUiQueryOutputs,
  pendingSteps
} from "../selectors/workflow_selectors";

export function useInputStateManagement(state: VulcanState, dispatch: Dispatch<VulcanAction>) {
  const {inputs, workflow, status, data, session} = state;

  useEffect(() => {
    if (workflow == null) return;
    const deletes: string[] = [];

    Object.keys(inputs).forEach(source => {
      if (!(source in workflow.dependencies_of_outputs)) {
        deletes.push(source);
        return;
      }

      const unfulfilled = findSourceDependencies(source, workflow).find(dependency => {
        return inputValueNonEmpty(inputs[dependency]);
      })

      if (unfulfilled) {
        deletes.push(source);
      }
    });

    console.log('deletes', deletes);

    if (deletes.length > 0) dispatch(removeInputs(deletes));
  }, [inputs, workflow]);

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

    console.log('newInputs', Object.keys(newInputs));

    if (Object.keys(newInputs).length > 0) dispatch(patchInputs(newInputs));
  }, [inputs, workflow, status, data, session]);
}