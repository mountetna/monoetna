import React from 'react';

import Icon from 'etna-js/components/icon';

const SaveCancel = ({onSave, onCancel}) => {
  return (
    <React.Fragment>
      <Icon className='approve' icon='save' onClick={onSave} />
      <Icon className='cancel' icon='ban' onClick={onCancel} />
    </React.Fragment>
  );
};

export default SaveCancel;
