import React, {useContext, useEffect, useState} from 'react';

import {VulcanContext} from '../../contexts/vulcan';
import {defaultInputValues} from '../../utils/workflow';
import {workflowByName} from '../../selectors/workflow';

import Link from 'etna-js/components/link';

import SessionManager from './session/session_manager';
import StepsList from './steps/steps_list';

export default function WorkflowManager({workflowName}) {
  const {
    workflows,
    setWorkflow,
    setPathIndex,
    setInputs,
    setSession,
    getLocalSession
  } = useContext(VulcanContext);

  useEffect(() => {
    if (
      workflows.workflows &&
      workflowByName({
        workflows: workflows.workflows,
        workflowName
      })
    ) {
      let selectedWorkflow = workflowByName({
        workflows: workflows.workflows,
        workflowName
      });
      setWorkflow(selectedWorkflow);

      getLocalSession(selectedWorkflow).then((session) => {
        if (null == session) {
          // Set the default input values
          setInputs(defaultInputValues(selectedWorkflow));
        } else {
          setSession(session);
        }
      });

      // first path is always the "work" path
      setPathIndex(0);
    }
  }, [workflows]);

  return (
    <div className='workflow-manager'>
      <SessionManager></SessionManager>
      <StepsList></StepsList>
    </div>
  );
}
