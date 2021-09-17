import React, {PropsWithChildren} from 'react';
import Icon from 'etna-js/components/icon';

export default function InputHelp({doc = '', children}: PropsWithChildren<{doc?: string}>) {
  return (
    <div className='input-help'>
      <div className='input-help-children-wrapper'>{children}</div>
      <div className='help-icon-wrapper' title={doc}>
        {doc ? (
          <Icon icon='question-circle' className='help-icon'/>
        ) : null}
      </div>
    </div>
  );
}
