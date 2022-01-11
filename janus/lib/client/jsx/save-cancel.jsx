import React from 'react';

import 'regenerator-runtime/runtime';

import Icon from 'etna-js/components/icon';
import useAsyncWork from 'etna-js/hooks/useAsyncWork';

const SaveCancel = ({onSave, onCancel}) => {
  const [_, safeOnSave] = useAsyncWork(onSave, {cancelWhenChange: []});
  return (
    <React.Fragment>
      <Icon className='approve' icon='save' onClick={safeOnSave} />
      <Icon className='cancel' icon='ban' onClick={onCancel} />
    </React.Fragment>
  );
};

export default SaveCancel;
