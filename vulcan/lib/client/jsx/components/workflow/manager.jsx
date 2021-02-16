import React, {useContext, useEffect} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from 'etna-js/actions/message_actions';

import {VulcanContext} from '../../contexts/vulcan';
import {getWorkflow} from '../../api/vulcan';
import Output from './output';
import Input from './input';

import StepsList from './steps/steps_list';

// Hardcode for now, since only one workflow
const WORKFLOW_NAME = 'umap';

export default function Manager() {
  const invoke = useActionInvoker();
  const {setWorkflow, setPathIndex} = useContext(VulcanContext);

  useEffect(() => {
    getWorkflow(WORKFLOW_NAME, true)
      .then((workflowDetails) => {
        setWorkflow(workflowDetails);

        // Hardcode the path for now, since only 1?
        setPathIndex(0);
      })
      .catch((e) => {
        invoke(showMessages(e));
      });
  }, []);

  return (
    <div className='workflow-manager'>
      <div className='workflow-header'>
        You are working on a {WORKFLOW_NAME.toUpperCase()} workflow:
      </div>
      <div className='step-wrapper'>
        <div className='step-nav-wrapper'>
          <StepsList></StepsList>
        </div>
        <div className='step-main-pane-wrapper'>
          <Input></Input>
          <Output></Output>
        </div>
      </div>
    </div>
  );
}
