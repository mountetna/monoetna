import React from 'react';
import {VulcanState} from "../reducers/vulcan_reducer";
import {
  createContext, Dispatch, MutableRefObject, PropsWithChildren, useCallback, useContext, useRef, useState
} from "react";
import {defaultSessionSyncHelpers} from "./session_sync";
import {useActionInvoker} from "etna-js/hooks/useActionInvoker";
import {showMessages} from "etna-js/actions/message_actions";
import {clearBufferedInput, setBufferedInput, setInputs, VulcanAction} from "../actions/vulcan_actions";
import {allSourcesForStepName} from "../selectors/workflow_selectors";
import {mapSome, Maybe} from "../selectors/maybe";
import {DataEnvelope} from "../components/workflow/user_interactions/inputs/input_types";

export const defaultInputStateManagement = {
  commitSessionInputChanges(stepName: string | null, inputs: DataEnvelope<Maybe<any>>) {
    return false;
  }
}

export const defaultBufferedInputs = {
  inputs: {} as DataEnvelope<Maybe<any>>,
  setInputs(inputs: DataEnvelope<Maybe<any>> | ((prev: DataEnvelope<Maybe<any>>) => DataEnvelope<Maybe<any>>)) {
  },
  commitInputs() {},
  cancelInputs() {},
}
export const BufferedInputsContext = createContext(defaultBufferedInputs);

export function WithBufferedInputs({
  children,
  commitSessionInputChanges,
  dispatch,
  stepName,
}: PropsWithChildren<{
  commitSessionInputChanges: typeof defaultInputStateManagement.commitSessionInputChanges,
  dispatch: Dispatch<VulcanAction>,
  stepName: string | null,
}>) {
  const inputsRef = useRef({} as DataEnvelope<Maybe<any>>);
  const [inputs, setInputsState] = useState(inputsRef.current);
  const hasInputs = Object.keys(inputs).length > 0;

  const setInputs = useCallback<typeof defaultBufferedInputs.setInputs>(inputs => {
    if (inputs instanceof Function) {
      inputsRef.current = inputs(inputsRef.current);
    } else {
      inputsRef.current = inputs;
    }

    if (Object.keys(inputsRef.current).length > 0){
      dispatch(setBufferedInput(stepName));
    } else {
      dispatch(clearBufferedInput(stepName));
    }
    setInputsState(inputsRef.current);
  }, [dispatch, stepName]);

  const cancelInputs = useCallback(() => {
    setInputs({});
  }, [setInputs])

  const commitInputs = useCallback(() => {
    if (commitSessionInputChanges(stepName, inputsRef.current)) {
      cancelInputs();
    }
  }, [cancelInputs, commitSessionInputChanges, stepName]);

  return <BufferedInputsContext.Provider value={{inputs, setInputs, commitInputs, cancelInputs}}>
    <div>
      {children}
    </div>
    { hasInputs ? <div>
      <button onClick={cancelInputs}>
        Cancel
      </button>
      <button onClick={commitInputs}>
        Commit
      </button>
    </div> : null }
  </BufferedInputsContext.Provider>
}

export function useInputStateManagement(invoke: ReturnType<typeof useActionInvoker>,
  dispatch: Dispatch<VulcanAction>,
  requestPoll: typeof defaultSessionSyncHelpers.requestPoll,
  stateRef: MutableRefObject<VulcanState>,
): typeof defaultInputStateManagement {
  const getErrors = useCallback((step: string | null) => stateRef.current.validationErrors
    .filter(([stepErr]) => step === stepErr)
    .map(([step, inputLabel, errors]) => {
      return errors.map((e: string) => `${inputLabel}: ${e}`);
    }).flat(), [stateRef]);

  const validateInputs = useCallback((step: string | null) => {
    const validationErrs = getErrors(step);
    if (validationErrs.length > 0) {
      invoke(showMessages(validationErrs));
    }

    return validationErrs.length === 0;
  }, [getErrors, invoke])

  const commitSessionInputChanges = useCallback<typeof defaultInputStateManagement.commitSessionInputChanges>((stepName, inputs) => {
    if (!validateInputs(stepName)) return false;
    const sources = allSourcesForStepName(stepName, stateRef.current.workflow);
    const newInputs = {...stateRef.current.session.inputs};
    sources.forEach(source => {
      mapSome(inputs[source] || null, inner => newInputs[source] = inner);
    })
    dispatch(setInputs(newInputs))
    requestPoll();
    return true;
  }, [dispatch, requestPoll, stateRef, validateInputs]);

  return {
    commitSessionInputChanges,
  };
}
