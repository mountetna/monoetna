import React, {useContext} from 'react';
import 'regenerator-runtime/runtime';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {pushLocation} from 'etna-js/actions/location_actions';

import {VulcanContext} from '../contexts/vulcan_context';
import Card from '../components/dashboard/card';

export default function Dashboard() {
  const invoke = useActionInvoker();
  let {state} = useContext(VulcanContext);
  const {workflows} = state;

  console.log({workflows});

  if (!workflows || workflows.length === 0) return null;

  function handleOnClick(workflow) {
    invoke(pushLocation(`/workflow/${workflow.name}`));
  }

  return (
    <main className='vulcan-dashboard'>
      {workflows.map((w, ind) => {
        return (
          <Card
            workflow={w}
            key={ind}
            onClick={() => {
              handleOnClick(w);
            }}
          />
        );
      })}
    </main>
  );
}
