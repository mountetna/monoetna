import React, {useState, useContext, useEffect} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from 'etna-js/actions/message_actions';
import Icon from 'etna-js/components/icon';
import Link from 'etna-js/components/link';

import {VulcanContext} from '../../../contexts/vulcan_context';
import InputFeed from './input_feed';
import OutputFeed from './output_feed';
import {allInputsDefined, workflowName} from '../../../utils/workflow';

export default function SessionManager() {
  const invoke = useActionInvoker();
  const context = useContext(VulcanContext);
  const {state} = context;

  // const [complete, setComplete] = useState(null);
  // const [firstRun, setFirstRun] = useState(true);

  // const complete = state.workflow

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
        <span className='session-workflow-name'>{workflowName(workflow)}</span>
        {workflow.vignette && (
          <Link link={ROUTES.workflow_vignette(workflowName(workflow))}>
            <Icon className='vignette' icon='book' />
          </Link>
        )}
        <Icon
          className='run'
          disabled={!complete}
          title='Run workflow'
          onClick={runWorkflow}
          icon='play'
        />
      </div>
      <div className='session-feed-container'>
        <InputFeed></InputFeed>
        <OutputFeed></OutputFeed>
      </div>
    </div>
  );
}
