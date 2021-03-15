import React, {useEffect, useContext} from 'react';
import 'regenerator-runtime/runtime';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {pushLocation} from 'etna-js/actions/location_actions';

import {VulcanContext} from '../contexts/vulcan';
import Card from '../components/dashboard/card';

export default function Dashboard() {
  const invoke = useActionInvoker();
  let {workflows} = useContext(VulcanContext);

  if (!workflows.workflows || workflows.workflows.length === 0) return null;

  function handleOnClick(workflow) {
    invoke(pushLocation(`/workflow/${workflow.name.replace('.cwl', '')}`));
  }

  return (
    <main className='vulcan-dashboard'>
      {workflows.workflows.map((w, ind) => {
        return (
          <Card
            workflow={w}
            key={ind}
            onClick={() => {
              handleOnClick(w);
            }}
          ></Card>
        );
      })}
    </main>
  );
}
