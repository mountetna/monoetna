import React from 'react';
import Icon from 'etna-js/components/icon';

export default function InputHelp({input, children}) {
  return (
    <div className='input-help'>
      {icon.doc ? (
        <Icon icon='question-circle' className='help-icon'></Icon>
      ) : null}
      {children}
    </div>
  );
}
