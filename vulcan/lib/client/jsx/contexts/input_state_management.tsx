import {VulcanState} from "../reducers/vulcan_reducer";
import {Dispatch, MutableRefObject, useCallback} from "react";
import {defaultSessionSyncHelpers} from "./session_sync";
import {useActionInvoker} from "etna-js/hooks/useActionInvoker";
import {showMessages} from "etna-js/actions/message_actions";
import {setBufferedInput, setInputs, VulcanAction} from "../actions/vulcan_actions";
import { allSourcesForStepName } from "../selectors/workflow_selectors";
import {mapSome} from "../selectors/maybe";
import {Workflow, WorkflowStep} from "../api_types";

export const defaultInputStateManagement = {
  commitSessionInputChanges(stepName: string | null) {
  },
  cancelInputChanges(stepName: string | null) {
  },
}

export function useInputStateManagement(
  invoke: ReturnType<typeof useActionInvoker>,
  dispatch: Dispatch<VulcanAction>,
  requestPoll: typeof defaultSessionSyncHelpers.requestPoll,
  stateRef: MutableRefObject<VulcanState>,
): typeof defaultInputStateManagement {
  const getErrors = useCallback(() => stateRef.current.validationErrors
    .map(([inputLabel, errors]) => {
      return errors.map((e: string) => `${inputLabel}: ${e}`);
    }).flat(), [stateRef]);

  const validateInputs = useCallback(() => {
    const validationErrs = getErrors();
    if (validationErrs.length > 0) {
      invoke(showMessages(validationErrs));
    }

    return validationErrs.length === 0;
  }, [getErrors, invoke])

  const cancelInputChanges = useCallback<typeof defaultInputStateManagement.cancelInputChanges>(stepName => {
    const sources = allSourcesForStepName(stepName, stateRef.current.workflow);
    const newBuffered = {...stateRef.current.bufferedInputValues};
    sources.forEach(source => {
      delete newBuffered[source];
    })

    dispatch(setBufferedInput(newBuffered));
  }, [dispatch, stateRef]);

  const commitSessionInputChanges = useCallback<typeof defaultInputStateManagement.commitSessionInputChanges>(stepName => {
    if (!validateInputs()) return;
    const sources = allSourcesForStepName(stepName, stateRef.current.workflow);
    const newInputs = {...stateRef.current.session.inputs};
    sources.forEach(source => {
      mapSome(stateRef.current.bufferedInputValues[source], inner => newInputs[source] = inner);
    })
    dispatch(setInputs(newInputs))
    requestPoll();
  }, [dispatch, requestPoll, stateRef, validateInputs]);

  return {
    cancelInputChanges,
    commitSessionInputChanges,
  };
}
