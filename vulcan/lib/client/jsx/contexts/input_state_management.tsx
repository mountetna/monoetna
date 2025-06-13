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
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {dismissMessages, showMessages} from 'etna-js/actions/message_actions';
import {
  clearBufferedInput,
  clearAutoPassStep,
  setBufferedInput,
  setUIValues,
  setAutoPassStep,
  VulcanAction
} from '../actions/vulcan_actions';
import {
  Maybe,
} from '../selectors/maybe';
import {VulcanContext} from './vulcan_context';

import FormControlLabel from '@material-ui/core/FormControlLabel';
import Grid from '@material-ui/core/Grid';
import Switch from '@material-ui/core/Switch';
import { Workspace } from '../api_types';
import { DataEnvelope } from '../components/workspace/ui_definitions/input_types';
import FlatButton from 'etna-js/components/flat-button';

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
  cancelValueUpdates() {},
  showError(e: string) {}
};
export const BufferedInputsContext = createContext(defaultBufferedInputs);

export function WithBufferedInputs({
  children,
  commitSessionInputChanges,
  dispatch,
  invoke,
  stepName
}: PropsWithChildren<{
  commitSessionInputChanges: typeof defaultInputStateManagement.commitSessionInputChanges;
  dispatch: Dispatch<VulcanAction>;
  invoke: ReturnType<typeof useActionInvoker>,
  stepName: string;
}>) {

  const {stateRef, state} = useContext(VulcanContext);
  const valuesRef = useRef({} as DataEnvelope<Maybe<any>>);
  const [values, setValuesState] = useState(valuesRef.current);
  // We must assume a newly rendered UI might be filled in automatically with a valid default, even if the user has not touched it. 
  const hasNewValues = Object.keys(values).length > 0 ||
    !((stepName in stateRef.current.status.params) || (stepName in stateRef.current.status.ui_contents));

  const cancelValueUpdates: any = useCallback(() => {
    // // eslint-disable-next-line
    setValues({});
    // // @ts-ignore
    // // eslint-disable-next-line
  }, [setValues]);

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

  const showError = useCallback((e: string) => {
    invoke(dismissMessages());
    invoke(showMessages([`Error in step ${stepName}: ${e}`]))
  }, [])

  const commit_reset_buttons = hasNewValues ? (
    <Grid container direction='row' className={'ui-controls'} justifyContent='flex-end'>
      <Grid item>
        <FlatButton
          className={'control-btn clear'}
          icon='eraser'
          label='Clear'
          title='Coming Soon, Reset to original (blank or defaulted) state'
          onClick={() => {}}
          disabled={true}
        />
      </Grid>
      <Grid item>
        <FlatButton
          className={'control-btn reset'}
          icon='backward'
          label='Back'
          title='Reset to last saved state'
          onClick={cancelValueUpdates}
        />
      </Grid>
      <Grid item>
        <FlatButton
          className={'control-btn confirm'}
          icon='save'
          label='Save Choices'
          title={state.attemptingToRun || state.isRunning ?
            'Cannot save new choices to workspace while work is running' :
            state.isSyncing ?
            'Awaiting sync from remote workspace' :
            'Save choices to workspace on computation server'}
          onClick={commitValueUpdates}
          disabled={state.isSyncing || state.isRunning}
        />
      </Grid>
    </Grid>
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
          disabled={state.isSyncing}
        />
      }
      label='Auto-Confirm in future'
      labelPlacement='start'
    />
  ) : null;

  const controls_below =
    isPassableUIStep(stepName, stateRef.current.workspace) || hasNewValues ? (
      <Grid container style={{width: 'auto'}} justifyContent='flex-end'>
        {autopass_switch}
        {commit_reset_buttons}
      </Grid>
    ) : null;

  return (
    <BufferedInputsContext.Provider
      value={{values, setValues, commitValueUpdates, cancelValueUpdates, showError}}
    >
      <div>{children}</div>
      {controls_below}
    </BufferedInputsContext.Provider>
  );
}

export function useInputStateManagement(
  invoke: ReturnType<typeof useActionInvoker>,
  dispatch: Dispatch<VulcanAction>,
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
      // Validate current input set, ending here if not valid
      if (!validateInputs(stepName)) return false;
      // Push to params/ui_contents in status, and mark step as needing to by sync'd to server by workspace-manager
      dispatch(setUIValues(inputs, stepName))
      return true;
    },
    [dispatch, validateInputs]
  );

  return {
    commitSessionInputChanges
  };
}
