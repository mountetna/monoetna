import React, {useMemo, useState, useCallback} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from 'etna-js/actions/message_actions';
import {addTemplatesAndDocuments} from 'etna-js/actions/magma_actions';
import {useModal} from 'etna-js/components/ModalDialogContainer';

const useMagmaActions = () => {
  const invoke = useActionInvoker();
  const {dismissModal} = useModal();

  const executeAction = useCallback(
    (action: any, closeModal = true) => {
      if (closeModal) dismissModal();
      return action
        .then(({models}: {models: any}) => {
          invoke(addTemplatesAndDocuments(models));
        })
        .catch((err: any) => {
          invoke(showMessages(err));
          throw err;
        });
    },
    [invoke, addTemplatesAndDocuments, showMessages]
  );
  return {
    executeAction
  };
};

export default useMagmaActions;
