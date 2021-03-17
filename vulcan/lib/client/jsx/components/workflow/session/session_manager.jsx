import React, {useState, useContext, useEffect} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from 'etna-js/actions/message_actions';
import Icon from 'etna-js/components/icon';
import Link from 'etna-js/components/link';

import {submit} from '../../../api/vulcan';
import {VulcanContext} from '../../../contexts/vulcan';
import {allInputsDefined, workflowName} from '../../../utils/workflow';
import SessionFeed from './session_feed';

import PrimaryInputs from './primary_inputs';

export default function SessionManager({name}) {
  const invoke = useActionInvoker();
  const context = useContext(VulcanContext);
  const {workflow, session, calculating, setCalculating} = context;
  const [complete, setComplete] = useState(null);
  const [firstRun, setFirstRun] = useState(true);

  useEffect(() => {
    if (
      workflow &&
      workflow.name &&
      allInputsDefined(workflow, session.inputs) &&
      firstRun
    ) {
      setCalculating(true);
      setFirstRun(false);

      // See how much work can be done when the page loads
      submit(context)
        .then(() => {
          setCalculating(false);
        })
        .catch((e) => {
          console.error(e);
          invoke(showMessages(e));
        });
    }
  }, [session.inputs]);

  useEffect(() => {
    if (
      workflow &&
      session &&
      session.inputs &&
      Object.keys(session.inputs).length > 0
    ) {
      setComplete(allInputsDefined(workflow, session.inputs));
    }
  }, [workflow, session]);

  function runWorkflow() {
    setCalculating(true);
    submit(context)
      .then(() => {
        setCalculating(false);
        setComplete(allInputsDefined(workflow, session.inputs));
      })
      .catch((e) => {
        console.error(e);
        invoke(showMessages(e));
      });
  }
  return (
    <div className='session-manager'>
      <div className='session-header'>
        <span className='session-workflow-name'>
        { workflowName(workflow) }
        </span>
        <Link link={ROUTES.workflow_vignette(workflowName(workflow))}>
          <Icon className='vignette' icon='book'/>
        </Link>
        <Icon
          className='run'
          disabled={!complete || calculating}
          title='Run workflow'
          onClick={runWorkflow}
          icon='play'/>
      </div>
      <div className='scroll-window'>
        <PrimaryInputs></PrimaryInputs>
        <SessionFeed></SessionFeed>
      </div>
    </div>
  );
}
