import React, {useState, useContext, useEffect} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from 'etna-js/actions/message_actions';

import {getSession, submit} from '../../../api/vulcan';
import {VulcanContext} from '../../../contexts/vulcan';
import {defaultInputValues} from '../../../selectors/workflow_selector';
import {STATUS} from '../../../models/steps';

export default function StartSession() {
  // Placeholder for when user can select a session
  //   or continue a past session.
  const invoke = useActionInvoker();
  const {workflow, pathIndex, setStepIndex, setSession, setStatus} = useContext(
    VulcanContext
  );

  function handleOnClick() {
    // Leave this as getSession (without default inputs), since
    //   we may provide a session key in the future to fetch
    //   a persisted session.
    getSession(workflow.name)
      .then((response) => {
        setSession(response.session);
        setStatus(response.status);

        // Submit any default inputs to kick off the calculations?
        // This should launch cells that have no inputs or
        //   that have all default inputs...
        return submit(
          workflow.name,
          defaultInputValues(workflow),
          response.session.key
        );
      })
      .then((response) => {
        // Let's update stepIndex with any uncompleted
        //   work, based on the updated status.
        setSession(response.session);
        setStatus(response.status);

        let stepIndex = response.status[pathIndex].findIndex((step) => {
          return STATUS.COMPLETE !== step.status;
        });

        setStepIndex(stepIndex);
      })
      .catch((e) => {
        console.error(e);
        invoke(showMessages([e]));
      });
  }
  return (
    <section className='start-session'>
      <button className='large' onClick={handleOnClick}>
        Start New Session
      </button>
    </section>
  );
}
