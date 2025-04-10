import React, {
  useContext,
  useState,
  useMemo
} from 'react';

import {VulcanContext} from '../../../contexts/vulcan_context';
import StepIcon from './step_elements/step_icon';
import {labelOfStepOrGroupedStep} from '../../../selectors/workflow_selectors';
import {statusOfStep} from '../../../selectors/workflow_selectors';
import {
  STATUS,
  WorkspaceStep
} from '../../../api_types';
import {
  WorkspaceStepGroup,
  bindInputSpecification,
  BoundInputSpecification,
  getInputSpecifications,
  InputSpecification
} from '../ui_definitions/input_types';

import Card from '@material-ui/core/Card';
import Typography from '@material-ui/core/Typography';
import Collapse from '@material-ui/core/Collapse';
import IconButton from '@material-ui/core/IconButton';
import Grid from '@material-ui/core/Grid';

import ExpandMoreIcon from '@material-ui/icons/ExpandMore';
import ExpandLessIcon from '@material-ui/icons/ExpandLess';

import UserInput from './user_input';
import {useWorkspace} from '../../../contexts/workspace_context';
import {
  BufferedInputsContext,
  WithBufferedInputs
} from '../../../contexts/input_state_management';
import { useInputFeedStyles } from '../input_feed';

function StepUIInner({
  specs,
  hideLabel
}: {
  specs: InputSpecification[];
  hideLabel: boolean;
}) {
  const {state} = useContext(VulcanContext);
  const {status} = state;
  const {workspace} = useWorkspace();
  const {values, setValues, showError} = useContext(BufferedInputsContext);

  const stepInputs: BoundInputSpecification[] = useMemo(
    () =>
      specs.map((spec) =>
        bindInputSpecification(
          spec,
          workspace.vulcan_config,
          status.last_params,
          status.file_contents,
          status.params,
          status.ui_contents,
          values,
          setValues,
          showError
        )
      ),
    [specs, workspace.vulcan_config, status.last_params, status.file_contents, status.params, status.ui_contents, values, setValues]
  );

  return (
    <React.Fragment>
      {stepInputs.map((input, index) => {
        return <UserInput input={input} hideLabel={hideLabel} key={index} />;
      })}
    </React.Fragment>
  );
}

function StepUIOuter({
  step,
  hideLabel = true
}: {
  step: WorkspaceStep;
  hideLabel: boolean;
}) {
  const {dispatch, commitSessionInputChanges, useActionInvoker} = useContext(VulcanContext);
  const {workspace} = useWorkspace();
  const specs = useMemo(
    () => getInputSpecifications(step, workspace),
    [step, workspace]
  );

  return (
    <React.Fragment>
      <WithBufferedInputs
        commitSessionInputChanges={commitSessionInputChanges}
        dispatch={dispatch}
        invoke={useActionInvoker()}
        stepName={step.name}
      >
        <StepUIInner specs={specs} hideLabel={hideLabel} />
      </WithBufferedInputs>
    </React.Fragment>
  );
}

const collator = new Intl.Collator(undefined, {
  numeric: true, sensitivity: 'base'
});

export default function GroupedStepUI({
  group
}: {
  group: WorkspaceStepGroup;
}) {
  const {state} = useContext(VulcanContext);

  const {status} = state;

  const allInnerStatus = useMemo(
    () => group.steps.map((step) => statusOfStep(step, status, state.workspace)),
    [group, status]
  );
  const allStepsComplete = useMemo(
    () => allInnerStatus.every((s) => s && s.status === STATUS.COMPLETE),
    [allInnerStatus]
  );
  const [open, setOpen] = useState(!allStepsComplete);

  const hasValidationErrors = group.steps.some((step) =>
    state.validationErrors.some(([stepName]) => stepName === step.name)
  );

  const classes = useInputFeedStyles();

  const label = labelOfStepOrGroupedStep(group);

  let stepInputs = useMemo(() => group.steps
    .sort((a, b) => collator.compare(a.label || a.name, b.label || b.name))
    .map(step => <StepUIOuter key={step.name} step={step} hideLabel={true}/>
  ), [group.steps]);
  
  return (
    <Card
      elevation={0}
      className={`${classes.card} ${
        hasValidationErrors ? classes.error : ''
      } step-user-input`}
    >
      <Grid
        className={classes.header}
        justifyContent='space-between'
        container
        onClick={() => setOpen(!open)}
      >
        <Grid item container style={{width: 'auto'}}>
          <StepIcon step={group} />
          <Typography className={classes.label}>{label}</Typography>
        </Grid>
        <IconButton size='small' onClick={() => setOpen(!open)}>
          {open ? (
            <ExpandLessIcon fontSize='small' />
          ) : (
            <ExpandMoreIcon fontSize='small' />
          )}
        </IconButton>
      </Grid>
      <Collapse className='step-user-input-inputs' in={open}>
        {stepInputs}
      </Collapse>
    </Card>
  );
}
