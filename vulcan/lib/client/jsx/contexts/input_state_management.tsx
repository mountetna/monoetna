import React, {useEffect, useMemo} from 'react';
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
import {defaultSessionSyncHelpers} from './session_sync_while_running';
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
import {stepOutputData, stepOutputMapping} from '../selectors/workflow_selectors';
import {
  mapSome,
  Maybe,
  some,
} from '../selectors/maybe';
import {VulcanContext} from './vulcan_context';

import Button from '@material-ui/core/Button';
import {FormControlLabel, Grid, Switch} from '@material-ui/core';
import { Workspace } from '../api_types';
import { DataEnvelope } from '../components/ui_components';

export const defaultInputStateManagement = {
  commitSessionInputChanges(
    stepName: string,
    inputs: DataEnvelope<Maybe<any>>
  ) {
    return false;
  }
};

export function isPassableUIStep(
  stepName: string,
  workspace: Workspace | null
) {
  if (workspace && Object.keys(workspace.vulcan_config).includes(stepName)) {
    const this_step = workspace?.vulcan_config[stepName];
    if (this_step?.doc?.startsWith('SKIPPABLE')) return true;
  }
  return false;
}

export const defaultBufferedInputs = {
  values: {} as DataEnvelope<Maybe<any>>,
  setValues(
    values:
      | DataEnvelope<Maybe<any>>
      | ((prev: DataEnvelope<Maybe<any>>) => DataEnvelope<Maybe<any>>)
  ) {},
  commitValueUpdates() {},
  cancelValueUpdates() {}
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
  stepName: string;
}>) {
  const {stateRef, state} = useContext(VulcanContext);
  const valuesRef = useRef({} as DataEnvelope<Maybe<any>>);
  const [values, setValuesState] = useState(valuesRef.current);
  const hasValues = Object.keys(values).length > 0;

  const cancelValueUpdates: any = useCallback(() => {
    if (stateRef.current.workspace!=null) {
      const config = (stateRef.current.workspace as Workspace).vulcan_config[stepName];
      const origValue = stepOutputData(
        stepName,
        stepOutputMapping(config),
        {},
        stateRef.current.status.params,
        stateRef.current.status.ui_contents,
        config['default'])

        // // eslint-disable-next-line
        setValues(origValue);
        // // @ts-ignore
        // // eslint-disable-next-line
    }
  }, [setValues, stepName, stateRef.current.workspace, stateRef.current.status.params, stateRef.current.status.ui_contents]);

  const setValues: any = useCallback<typeof defaultBufferedInputs.setValues>(
    (values) => {
      if (values instanceof Function) {
        valuesRef.current = values(valuesRef.current);
      } else {
        valuesRef.current = values;
      }

      if (Object.keys(valuesRef.current).length > 0) {
        if (!stateRef.current.bufferedSteps.includes(stepName)) {
          dispatch(setBufferedInput(stepName));
        }
        // // Check / Initiate auto-pass attempt stepUI.
        // if (
        //   stepName != null &&
        //   stateRef.current.autoPassSteps.includes(stepName) &&
        //   Object.keys(stateRef.current.status.ui_contents).filter((val) =>
        //     val.includes(stepName)
        //   ).length < 1
        // ) {
        //   if (commitSessionInputChanges(stepName, valuesRef.current)) {
        //     dispatch(setRunTrigger(stepName));
        //   }
        // }
      } else {
        if (stateRef.current.bufferedSteps.includes(stepName)) {
          dispatch(clearBufferedInput(stepName));
        }
      }
      setValuesState(valuesRef.current);
    },
    [dispatch, stateRef, stepName, cancelValueUpdates, commitSessionInputChanges]
  );

  // useEffect(() => {
  //   // Initiate auto-pass attempt, primary inputs.
  //   if (
  //     stepName === null &&
  //     stateRef.current.autoPassSteps.includes(stepName)
  //   ) {
  //     if (commitSessionInputChanges(stepName, valuesRef.current)) {
  //       cancelValueUpdates();
  //       dispatch(setRunTrigger(stepName));
  //     }
  //     dispatch(clearAutoPassStep(null));
  //   }
  // }, [stateRef.current.autoPassSteps, cancelValueUpdates, commitSessionInputChanges]);

  const commitValueUpdates = useCallback(() => {
    if (commitSessionInputChanges(stepName, valuesRef.current)) {
      cancelValueUpdates();
    }
  }, [cancelValueUpdates, commitSessionInputChanges, stepName]);

  function setAutoPass(event: any, checked: boolean) {
    if (checked) {
      dispatch(setAutoPassStep(stepName));
    } else {
      dispatch(clearAutoPassStep(stepName));
    }
  }

  const commit_reset_buttons = hasValues ? (
    <div className='reset-or-commit-inputs'>
      <Button onClick={cancelValueUpdates} disabled={!!state.pollingState}>
        Reset
      </Button>
      <Button
        onClick={commitValueUpdates}
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
    isPassableUIStep(stepName, stateRef.current.workspace) || hasValues ? (
      <Grid container style={{width: 'auto'}} justifyContent='flex-end'>
        {autopass_switch}
        {commit_reset_buttons}
      </Grid>
    ) : null;

  return (
    <BufferedInputsContext.Provider
      value={{values, setValues, commitValueUpdates, cancelValueUpdates}}
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
      // Set ui_contents in status
      dispatch(setUIValues(inputs, stepName))
      // Validate current input set, ending here if not valid
      if (!validateInputs(stepName)) return false;
      // Push for this step only
      requestPoll(true, false, stepName);
      return true;
    },
    [dispatch, requestPoll, validateInputs]
  );

  return {
    commitSessionInputChanges
  };
}
