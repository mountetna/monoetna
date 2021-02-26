import React, {useState, useContext, useEffect} from 'react';
import * as _ from 'lodash';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from 'etna-js/actions/message_actions';
import Icon from 'etna-js/components/icon';

import {getSession, submit} from '../../../api/vulcan';
import {VulcanContext} from '../../../contexts/vulcan';
import {allInputsDefined, defaultInputValues} from '../../../utils/workflow';
import SessionFeed from './session_feed';

import PrimaryInputs from './primary_inputs';

export default function SessionManager() {
  // Placeholder for when user can select a session
  //   or continue a past session.
  const invoke = useActionInvoker();
  const context = useContext(VulcanContext);
  const {
    workflow,
    session,
    setSession,
    setStatus,
    setInputs,
    setCalculating
  } = context;
  const [complete, setComplete] = useState(null);
  const [firstRun, setFirstRun] = useState(true);

  useEffect(() => {
    // Leave this as getSession (without default inputs), since
    //   we may provide a session key in the future to fetch
    //   a persisted session.
    if (workflow && workflow.name) {
      getSession(workflow.name)
        .then((response) => {
          setSession(response.session);
          setStatus(response.status);

          // Set the default input values in the session
          setInputs(defaultInputValues(workflow));
        })
        .catch((e) => {
          console.error(e);
          invoke(showMessages(e));
        });
    }
  }, [workflow]);

  useEffect(() => {
    if (workflow && session && session.inputs) {
      setComplete(allInputsDefined(workflow, session.inputs));
    }
  }, [workflow, session]);

  function handleOnClick() {
    setCalculating(true);
    setFirstRun(false);
    submit(context)
      .then(() => {
        setCalculating(false);
      })
      .catch((e) => {
        console.error(e);
        invoke(showMessages([e]));
      });
  }
  return (
    <div className='session-manager'>
      <div className='start-btn-container'>
        <button
          disabled={!complete}
          className='start-button'
          onClick={handleOnClick}
        >
          Run
          <Icon icon='play' className='small'></Icon>
        </button>
      </div>
      <div className='scroll-window'>
        <PrimaryInputs></PrimaryInputs>
        {!firstRun ? <SessionFeed></SessionFeed> : null}
      </div>
    </div>
  );
}
