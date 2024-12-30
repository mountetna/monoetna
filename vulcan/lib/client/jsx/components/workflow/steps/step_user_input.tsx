import React, {useContext, useMemo} from 'react';
import {VulcanContext} from '../../../contexts/vulcan_context';
import UserInput from '../user_interactions/inputs/user_input';
import {WorkspaceStep} from '../../../api_types';
import {
  bindInputSpecification,
  BoundInputSpecification,
  getInputSpecifications,
  InputSpecification
} from '../user_interactions/inputs/input_types';
import {useWorkspace} from '../../../contexts/workspace_context';
import {
  BufferedInputsContext,
  WithBufferedInputs
} from '../../../contexts/input_state_management';

function StepUserInputInner({
  specs,
  hideLabel
}: {
  specs: InputSpecification[];
  hideLabel: boolean;
}) {
  const {state} = useContext(VulcanContext);
  const {status} = state;
  const {workspace} = useWorkspace();
  const {values, setValues} = useContext(BufferedInputsContext);

  const stepInputs: BoundInputSpecification[] = useMemo(
    () =>
      specs.map((spec) =>
        bindInputSpecification(
          spec,
          workspace.steps,
          workspace.vulcan_config,
          status.last_params,
          status.file_contents,
          status.params,
          status.ui_contents,
          values,
          setValues
        )
      ),
    [specs, workspace, status, values, setValues]
  );

  return (
    <React.Fragment>
      {stepInputs.map((input, index) => {
        return <UserInput input={input} hideLabel={hideLabel} key={index} />;
      })}
    </React.Fragment>
  );
}

export default function StepUserInput({
  step,
  hideLabel = true
}: {
  step: WorkspaceStep;
  hideLabel: boolean;
}) {
  const {dispatch, commitSessionInputChanges} = useContext(VulcanContext);
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
        stepName={step.name}
      >
        <StepUserInputInner specs={specs} hideLabel={hideLabel} />
      </WithBufferedInputs>
    </React.Fragment>
  );
}
