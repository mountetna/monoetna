import React from 'react';
import Icon from 'etna-js/components/icon';

export default function InputHelp({input, children}) {
  return (
    <div className='input-help'>
      <div className='input-help-children-wrapper'>{children}</div>
      <div className='help-icon-wrapper' title={input.doc}>
        {input.doc ? (
          <Icon icon='question-circle' className='help-icon'></Icon>
        ) : null}
      </div>
    </div>
  );
}
