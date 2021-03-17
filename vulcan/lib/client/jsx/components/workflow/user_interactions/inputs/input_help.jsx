import React from 'react';
import Icon from 'etna-js/components/icon';

export default function InputHelp({input, children}) {
  return (
    <div className='input-help'>
      <div>{children}</div>
      {input.doc ? (
        <div className='help-icon-wrapper' title={input.doc}>
          <Icon icon='question-circle' className='help-icon'></Icon>
        </div>
      ) : null}
    </div>
  );
}
