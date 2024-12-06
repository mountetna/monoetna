import React, {useEffect} from 'react';
import {VulcanState} from '../reducers/vulcan_reducer';
import {
  createContext,
  Dispatch,
  MutableRefObject,
  PropsWithChildren,
  useCallback,
  useContext,
  useRef,
  useState
} from 'react';
import {defaultSessionSyncHelpers} from './session_sync';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {dismissMessages, showMessages} from 'etna-js/actions/message_actions';
import {
  clearBufferedInput,
  clearAutoPassStep,
  setBufferedInput,
  setUIValues,
  setRunTrigger,
  setAutoPassStep,
  VulcanAction
} from '../actions/vulcan_actions';
import {allSourcesForStepName} from '../selectors/workflow_selectors';
import {
  mapSome,
  Maybe,
  some,
} from '../selectors/maybe';
import {DataEnvelope} from '../components/workflow/user_interactions/inputs/input_types';
import {VulcanContext} from './vulcan_context';

import Button from '@material-ui/core/Button';
import {FormControlLabel, Grid, Switch} from '@material-ui/core';
import { Workspace } from '../api_types';

export const defaultInputStateManagement = {
  commitSessionInputChanges(
    stepName: string | null,
    inputs: DataEnvelope<Maybe<any>>
  ) {
    return false;
  }
};

export function isPassableUIStep(
  stepName: string | null,
  workspace: Workspace | null
) {
  if (workspace && stepName != null && Object.keys(workspace.vulcan_config).includes(stepName)) {
    const this_step = workspace?.vulcan_config[stepName];
    if (this_step?.doc?.startsWith('SKIPPABLE')) return true;
  }
  return false;
}

export const defaultBufferedInputs = {
  inputs: {} as DataEnvelope<Maybe<any>>,
  setInputs(
    inputs:
      | DataEnvelope<Maybe<any>>
      | ((prev: DataEnvelope<Maybe<any>>) => DataEnvelope<Maybe<any>>)
  ) {},
  commitInputs() {},
  cancelInputs() {}
};
export const BufferedInputsContext = createContext(defaultBufferedInputs);

export function WithBufferedInputs({
  children,
  commitSessionInputChanges,
  dispatch,
  stepName
}: PropsWithChildren<{
  commitSessionInputChanges: typeof defaultInputStateManagement.commitSessionInputChanges;
  dispatch: Dispatch<VulcanAction>;
  stepName: string | null;
}>) {
  const {stateRef, state} = useContext(VulcanContext);
  const inputsRef = useRef({} as DataEnvelope<Maybe<any>>);
  const [inputs, setInputsState] = useState(inputsRef.current);
  const hasInputs = Object.keys(inputs).length > 0;

  const cancelInputs: any = useCallback(() => {
    // eslint-disable-next-line
    setInputs({});
    // @ts-ignore
    // eslint-disable-next-line
  }, [setInputs]);

  const setInputs: any = useCallback<typeof defaultBufferedInputs.setInputs>(
    (inputs) => {
      if (inputs instanceof Function) {
        inputsRef.current = inputs(inputsRef.current);
      } else {
        inputsRef.current = inputs;
      }

      if (Object.keys(inputsRef.current).length > 0) {
        if (!stateRef.current.bufferedSteps.includes(stepName)) {
          dispatch(setBufferedInput(stepName));
        }
        // Check / Initiate auto-pass attempt stepUI.
        if (
          stepName != null &&
          stateRef.current.autoPassSteps.includes(stepName) &&
          Object.keys(stateRef.current.status.ui_contents).filter((val) =>
            val.includes(stepName)
          ).length < 1
        ) {
          if (commitSessionInputChanges(stepName, inputsRef.current)) {
            cancelInputs();
            dispatch(setRunTrigger(stepName));
          }
        }
      } else {
        if (stateRef.current.bufferedSteps.includes(stepName)) {
          dispatch(clearBufferedInput(stepName));
        }
      }
      setInputsState(inputsRef.current);
    },
    [dispatch, stateRef, stepName, cancelInputs, commitSessionInputChanges]
  );

  useEffect(() => {
    if (
      !state.bufferedSteps.includes(stepName) &&
      Object.keys(inputs).length > 0
    ) {
      setInputs({});
    }
  }, [setInputs, state.bufferedSteps, inputs, stepName]);

  useEffect(() => {
    // Initiate auto-pass attempt, primary inputs.
    if (
      stepName === null &&
      stateRef.current.autoPassSteps.includes(stepName)
    ) {
      if (commitSessionInputChanges(stepName, inputsRef.current)) {
        cancelInputs();
        dispatch(setRunTrigger(stepName));
      }
      dispatch(clearAutoPassStep(null));
    }
  }, [stateRef.current.autoPassSteps, cancelInputs, commitSessionInputChanges]);

  const commitInputs = useCallback(() => {
    if (commitSessionInputChanges(stepName, inputsRef.current)) {
      cancelInputs();
    }
  }, [cancelInputs, commitSessionInputChanges, stepName]);

  function setAutoPass(event: any, checked: boolean) {
    if (checked) {
      dispatch(setAutoPassStep(stepName));
    } else {
      dispatch(clearAutoPassStep(stepName));
    }
  }

  const commit_reset_buttons = hasInputs ? (
    <div className='reset-or-commit-inputs'>
      <Button onClick={cancelInputs} disabled={!!state.pollingState}>
        Reset
      </Button>
      <Button
        onClick={commitInputs}
        style={{
          background: 'linear-gradient(135deg, #6e8efb, #a777e3)',
          textAlign: 'center',
          fontWeight: 'bold',
          fontSize: 14,
          border: '3px solid black'
          //outline: '#4CAF50 solid 2px'
        }}
        disabled={!!state.pollingState}
      >
        Confirm
      </Button>
    </div>
  ) : null;

  const autopass_switch = isPassableUIStep(
    stepName,
    stateRef.current.workspace
  ) ? (
    <FormControlLabel
      control={
        <Switch
          checked={stateRef.current.autoPassSteps.includes(stepName)}
          onChange={setAutoPass}
          color='primary'
          disabled={!!state.pollingState}
        />
      }
      label='Auto-Confirm in future'
      labelPlacement='start'
    />
  ) : null;

  const controls_below =
    isPassableUIStep(stepName, stateRef.current.workspace) || hasInputs ? (
      <Grid container style={{width: 'auto'}} justifyContent='flex-end'>
        {autopass_switch}
        {commit_reset_buttons}
      </Grid>
    ) : null;

  return (
    <BufferedInputsContext.Provider
      value={{inputs, setInputs, commitInputs, cancelInputs}}
    >
      <div>{children}</div>
      {controls_below}
    </BufferedInputsContext.Provider>
  );
}

export function useInputStateManagement(
  invoke: ReturnType<typeof useActionInvoker>,
  dispatch: Dispatch<VulcanAction>,
  requestPoll: typeof defaultSessionSyncHelpers.requestPoll,
  stateRef: MutableRefObject<VulcanState>
): typeof defaultInputStateManagement {
  const getErrors = useCallback(
    (step: string | null) =>
      stateRef.current.validationErrors
        .filter(([stepErr]) => step === stepErr)
        .map(([step, inputLabel, errors]) => {
          return errors.map((e: string) => `${inputLabel}: ${e}`);
        })
        .flat(),
    [stateRef]
  );

  const validateInputs = useCallback(
    (step: string | null) => {
      invoke(dismissMessages());
      const validationErrs = getErrors(step);
      if (validationErrs.length > 0) {
        invoke(showMessages(validationErrs));
      }

      return validationErrs.length === 0;
    },
    [getErrors, invoke]
  );

  const commitSessionInputChanges = useCallback<
    typeof defaultInputStateManagement.commitSessionInputChanges
  >(
    (stepName, inputs) => {
      if (!validateInputs(stepName)) return false;
      const sources = allSourcesForStepName(
        stepName,
        stateRef.current.workspace
      );
      const newInputs = {...stateRef.current.status.config_contents, ...stateRef.current.status.ui_contents};
      sources.forEach((source) => {
        if (source.indexOf('/') > -1) {
          const [source_step, source_name] = source.split('/');
          mapSome(inputs[source_name] || null, (inner) => (newInputs[source_step][source_name] = inner));
        } else {
          mapSome(inputs[source] || null, (inner) => (newInputs[source] = inner));
        }
      });
      dispatch(setUIValues(newInputs));
      requestPoll(false, some(stepName as string));
      return true;
    },
    [dispatch, requestPoll, stateRef, validateInputs]
  );

  return {
    commitSessionInputChanges
  };
}
