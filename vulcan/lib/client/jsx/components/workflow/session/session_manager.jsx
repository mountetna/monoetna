import React, {useState, useContext, useEffect} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from 'etna-js/actions/message_actions';
import Icon from 'etna-js/components/icon';

import {submit} from '../../../api/vulcan';
import {VulcanContext} from '../../../contexts/vulcan';
import {allInputsDefined} from '../../../utils/workflow';
import SessionFeed from './session_feed';

import PrimaryInputs from './primary_inputs';

export default function SessionManager() {
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

      // See how much work can be done when the page loads
      submit(context)
        .then(() => {
          setCalculating(false);
          setFirstRun(false);
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

  function handleOnClick() {
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
      <div className='start-btn-container'>
        <button
          disabled={!complete || calculating}
          className='start-button'
          onClick={handleOnClick}
        >
          Run
          <Icon icon='play' className='small'></Icon>
        </button>
      </div>
      <div className='scroll-window'>
        <PrimaryInputs></PrimaryInputs>
        <SessionFeed></SessionFeed>
      </div>
    </div>
  );
}
