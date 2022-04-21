import React, {useEffect} from 'react';
import {VulcanState} from "../reducers/vulcan_reducer";
import {
  createContext, Dispatch, MutableRefObject, PropsWithChildren, useCallback, useContext, useRef, useState
} from "react";
import {defaultSessionSyncHelpers} from "./session_sync";
import {useActionInvoker} from "etna-js/hooks/useActionInvoker";
import {dismissMessages, showMessages} from "etna-js/actions/message_actions";
import {clearBufferedInput, clearCommitTrigger, clearAutoPassStep, setBufferedInput, setInputs, setRunTrigger, setAutoPassStep, VulcanAction} from "../actions/vulcan_actions";
import {allSourcesForStepName} from "../selectors/workflow_selectors";
import {mapSome, Maybe, maybeOfNullable, some, withDefault} from "../selectors/maybe";
import {DataEnvelope} from "../components/workflow/user_interactions/inputs/input_types";
import {VulcanContext} from "./vulcan_context";

import Button from '@material-ui/core/Button';
import { FormControlLabel, Grid, Switch } from '@material-ui/core';

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
  const {stateRef, state} = useContext(VulcanContext);
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
      if (!stateRef.current.bufferedSteps.includes(stepName))
        dispatch(setBufferedInput(stepName));
    } else {
      if (stateRef.current.bufferedSteps.includes(stepName))
        dispatch(clearBufferedInput(stepName));
    }
    setInputsState(inputsRef.current);
  }, [dispatch, stateRef, stepName]);

  useEffect(() => {
    if (!state.bufferedSteps.includes(stepName) && Object.keys(inputs).length > 0) {
      setInputs({});
    }
  }, [setInputs, state.bufferedSteps, inputs, stepName]);

  const cancelInputs = useCallback(() => {
    setInputs({});
  }, [setInputs])

  const commitInputs = useCallback(() => {
    if (commitSessionInputChanges(stepName, inputsRef.current)) {
      cancelInputs();
    }
  }, [cancelInputs, commitSessionInputChanges, stepName]);

  function setAutoPass(event: any, checked: boolean) {
    if (checked) {
      dispatch(setAutoPassStep(stepName))
    } else {
      dispatch(clearAutoPassStep(stepName))
    }
  }

  function isPassable(stepName: string | null) {
    if (stepName!=null) {
      // ToDo: better way?
      const this_step = state.workflow?.steps[0].filter((val) => val.name == stepName)[0]
      if (this_step && this_step.doc!=null && this_step.doc.startsWith('SKIPPABLE')) return true
    }
    return false
  }

  // Catch auto-pass 'Commit' trigger
  useEffect(() => {
    if (state.triggerCommit.length>0 && state.triggerCommit[0] == stepName) {
      // Try to commit, add 'Run' trigger if works
      if (commitSessionInputChanges(stepName, inputsRef.current)) {
        cancelInputs();
        dispatch(setRunTrigger(stepName))
      }
      dispatch(clearCommitTrigger(stepName))
    }
  },[state.triggerCommit])

  const commit_rest_buttons = hasInputs ? <div>
    <Button onClick={cancelInputs} disabled={!!state.pollingState}>
      Reset
    </Button>
    <Button
      onClick={commitInputs}
      style={{
        background: "linear-gradient(135deg, #6e8efb, #a777e3)",
        textAlign: 'center',
        fontWeight: 'bold',
        fontSize: 14,
        border:'3px solid black',
        //outline: '#4CAF50 solid 2px'
        }}
      disabled={!!state.pollingState}>
      Commit
    </Button>
  </div> : null

  const autopass_switch = isPassable(stepName) ? <FormControlLabel
    className='auto-pass-inputs'
    control={
      <Switch
        checked={state.autoPassSteps.includes(stepName)}
        onChange={setAutoPass}
        color='primary'
        disabled={!!state.pollingState}
      />
    }
    label="Auto-Commit in future"
    labelPlacement='start'
  /> : null

  const controls_below = (isPassable(stepName) || hasInputs) ? (
    <Grid container style={{width: 'auto'}} justify="flex-end" className='reset-or-commit-inputs'>
      {autopass_switch}
      {commit_rest_buttons}
    </Grid>
   ) : null

  return <BufferedInputsContext.Provider value={{inputs, setInputs, commitInputs, cancelInputs}}>
    <div>
      {children}
    </div>
    { controls_below } 
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
    invoke(dismissMessages());
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
    requestPoll(false, some(stepName as string))
    return true;
  }, [dispatch, requestPoll, stateRef, validateInputs]);

  return {
    commitSessionInputChanges,
  };
}
