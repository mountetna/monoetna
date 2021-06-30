import React, {useState, useContext, useEffect, useCallback, useMemo, Dispatch} from 'react';

import ExpandMoreIcon from '@material-ui/icons/ExpandMore';
import ExpandLessIcon from '@material-ui/icons/ExpandLess';

import {VulcanContext} from '../../../contexts/vulcan_context';

import UserInput from '../user_interactions/inputs/user_input';
import {
  completedSteps, inputGroupName, isPendingUiQuery, pendingSteps, uiQueryOfStep, sortInputsByLabel, useMemoized
} from '../../../selectors/workflow_selectors';
import {InputSpecification} from "../user_interactions/inputs/input_types";
import {useWorkflow} from "../../../contexts/workflow_context";
import {Maybe} from "../../../selectors/maybe";
import {Input} from "@material-ui/core";

interface Props {
  inputs: InputSpecification[],
  groupName: string,
}

export default function InputGroup({inputs, groupName}: Props) {
  const {state} = useContext(VulcanContext);
  const {workflow} = useWorkflow();
  let {session, status, data} = state;
  const [open, setOpen] = useState(true);

  // let groupName = inputGroupName(inputs[0].name);

  const toggleInputs = useCallback(() => setOpen(!open), [setOpen, open]);

  useEffect(() => {
    if (inputs.length === 0) return;

    if (session && session.inputs) {
      // Automatically collapse the input group if there are pending or completed UI input steps.
      let completed = completedSteps(workflow, status).filter(({step}) => !!uiQueryOfStep(step));
      let nextUiSteps = pendingSteps(workflow, status)
        .filter(({step}) => isPendingUiQuery(step, status, data, session));

      if (completed.length > 0 || nextUiSteps.length > 0) setOpen(false);
    }
  }, [status, data, session, workflow, inputs.length]);

  const sortedInputs = useMemoized(sortInputsByLabel, inputs);

  return (<div className='inputs-pane'>
      <div className='header-wrapper'>
        <div onClick={toggleInputs} className={`inputs-pane-header toggle ${open ? 'open' : 'closed'}`}>
          <div className='title'>{groupName}</div>
          {open ? <ExpandLessIcon/> : <ExpandMoreIcon/>}
        </div>
        <div className='filler'/>
      </div>
      <div
        className={`primary-inputs-container items sliding-panel vertical ${open ? 'open' : 'closed'}`}
      >
        {sortedInputs.map((input, index) => {
          return (<UserInput
              input={input}
              key={index}
            />);
        })}
      </div>
    </div>);
}
