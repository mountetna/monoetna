import React, {useState, useContext, useEffect} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from 'etna-js/actions/message_actions';
import Icon from 'etna-js/components/icon';

import {getSession, submit} from '../../../api/vulcan';
import {VulcanContext} from '../../../contexts/vulcan';
import {
  allInputsDefined,
  defaultInputValues
} from '../../../selectors/workflow_selector';

import PrimaryInputs from './primary_inputs';

export default function SessionManager({name}) {
  // Placeholder for when user can select a session
  //   or continue a past session.
  const invoke = useActionInvoker();
  const {
    workflow,
    session,
    pathIndex,
    setStepIndex,
    setSession,
    setStatus,
    setInputs
  } = useContext(VulcanContext);
  const [complete, setComplete] = useState(null);

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
          invoke(showMessages([e]));
        });
    }
  }, [workflow]);

  useEffect(() => {
    if (workflow && session && session.inputs) {
      setComplete(allInputsDefined(workflow, session.inputs));
    }
  }, [workflow, session]);

  function handleOnClick() {
    submit(workflow.name, session.inputs, session.key)
      .then((response) => {
        setSession(response.session);
        setStatus(response.status);
      })
      .catch((e) => {
        console.error(e);
        invoke(showMessages([e]));
      });
  }
  return (
    <div className='start-session'>
      <div className='session-header'>
        <div className='session-title'>{name.toUpperCase()}</div>
        <div className='session-run'>
          { complete ? <Icon icon='play' title='Run' className='run' onClick={handleOnClick}/> : <span className='waiting'>Awaiting input...</span> }
        </div>
      </div>
      <PrimaryInputs></PrimaryInputs>
    </div>
  );
}
