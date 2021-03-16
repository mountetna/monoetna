// Framework libraries.
import React, {useContext, useEffect} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from 'etna-js/actions/message_actions';

import {getWorkflows} from '../api/vulcan';
import {VulcanContext} from '../contexts/vulcan';

export default function ContextManager({params, children}) {
  const invoke = useActionInvoker();
  let {setWorkflows, setCalculating} = useContext(VulcanContext);

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

  return <React.Fragment>{children}</React.Fragment>;
}
