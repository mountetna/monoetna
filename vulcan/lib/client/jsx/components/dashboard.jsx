import React, {useEffect, useContext} from 'react';
import 'regenerator-runtime/runtime';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from 'etna-js/actions/message_actions';

import {VulcanContext} from '../contexts/vulcan';
import {getWorkflows} from '../api/vulcan';

export default function Dashboard() {
  const invoke = useActionInvoker();
  let {workflows, setWorkflows, setCalculating} = useContext(VulcanContext);

  useEffect(() => {
    setCalculating(true);
    getWorkflows()
      .then((response) => {
        setWorkflows(response);
        setCalculating(false);
      })
      .catch((e) => {
        console.error(e);
        invoke(showMessages([e]));
      });
  }, []);

  if (!workflows.workflows || workflows.workflows.length === 0) return null;

  return (
    <main className='vulcan-dashboard'>
      {workflows.workflows.map((w) => {
        return <div>{w.name}</div>;
      })}
    </main>
  );
}
