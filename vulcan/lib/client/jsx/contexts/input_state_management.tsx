import {VulcanState} from "../reducers/vulcan_reducer";
import {Dispatch, MutableRefObject, useCallback, useEffect} from "react";
import {defaultSessionSyncHelpers} from "./session_sync";
import {useActionInvoker} from "etna-js/hooks/useActionInvoker";
import {showMessages} from "etna-js/actions/message_actions";
import {clearBufferedInput, removeInputs, setBufferedInput, setInputs, VulcanAction} from "../actions/vulcan_actions";
import {inputSourcesOfStep, stepOfStatus} from "../selectors/workflow_selectors";
import {Cancellable} from "etna-js/utils/cancellable";
import {defaultConfirmationHelpers} from "./confirmation";

export const defaultInputStateManagement = {
  commitSessionInputChanges() {
    return Promise.resolve(true);
  },

  startInputChange(source: string, context: Cancellable) {
    return Promise.resolve(true);
  },

  cancelInputChange() {
  },

  onInputChange(source: string, value: any) {
  }
}

export function useInputStateManagement(
  state: MutableRefObject<VulcanState>,
  invoke: ReturnType<typeof useActionInvoker>,
  dispatch: Dispatch<VulcanAction>,
  requestPoll: typeof defaultSessionSyncHelpers.requestPoll,
  confirm: typeof defaultConfirmationHelpers.confirm,
): typeof defaultInputStateManagement {
  const commitSessionInputChanges = useCallback(() => {
    const {curEditingInput, bufferedInputValue, session} = state.current;
    const {inputs} = session;

    if (curEditingInput == null || bufferedInputValue == null) {
      return Promise.resolve(false);
    }

    if (Object.keys(state.current.validationErrors).length > 0) {
      invoke(showMessages(Object.entries(state.current.validationErrors)
        .map(([inputName, validation]: [string, any]) => {
          let {
            inputLabel, errors
          }: { inputLabel: string; errors: string[] } = validation;
          return errors.map((e: string) => `${inputLabel}: ${e}`);
        })
        .flat()));

      return Promise.resolve(false);
    }

    // This should update the state.current synchronously as well.
    dispatch(setInputs({...inputs, [curEditingInput]: bufferedInputValue[0]}));

    const preCommitStepStatuses = state.current.status[0];
    const [refreshWithInputs] = requestPoll();

    return refreshWithInputs.then(({result, cancelled}) => {
      if (cancelled || !result) return false;
      const workflow = state.current.workflow;
      if (!workflow) return false;

      // Aggregate together the status hashes, then remove each step whose hash remains unchanged.
      // The remaining keys will be stale step names, which can be expanded to all inputs to drop.
      const statusHashes: {[k: string]: string} = {};
      preCommitStepStatuses.forEach(({name, hash}) => {
        statusHashes[name] = hash;
      })
      result.status[0].forEach(({name, hash}) => {
        if (name in statusHashes && statusHashes[name] === hash) delete statusHashes[name];
      })

      // Clear any steps that either no longer report in status, or whose hashes have changed.
      dispatch(removeInputs(Object.keys(statusHashes).reduce((acc, stepName) => {
        const step = stepOfStatus(stepName, workflow);
        if (step) acc.push(...inputSourcesOfStep(step));
        return acc;
      }, [] as string[])))

      return true;
    });
  }, [dispatch, invoke, requestPoll, state]);

  const startInputChange = useCallback((source: string, context: Cancellable) => {
    if (state.current.curEditingInput != null) {
      return confirm("You have uncommitted changes to another input, but you may only change one input at a time.  Proceed to reset that input, or cancel these changes?", context).then(confirmed => {
        if (confirmed) {
          dispatch(setBufferedInput(source, null));
          return true;
        }

        return false;
      });
    }

    dispatch(setBufferedInput(source, null));
    return Promise.resolve(true);
  }, [confirm, state, dispatch]);

  const cancelInputChange = useCallback(() => {
    dispatch(clearBufferedInput);
  }, [dispatch]);

  const onInputChange = useCallback((source: string, value: any) => {
    if (state.current.curEditingInput !== source) {
      throw new Error(`Concurrent edit to ${source}, a call to startInputChange is necessary.`)
    }

    setBufferedInput(source, [value]);
  }, [state]);

  return {
    commitSessionInputChanges,
    startInputChange,
    cancelInputChange,
    onInputChange,
  }
}
