import React from 'react';
import Icon from 'etna-js/components/icon';
import {InputSpecification} from "./types";

export default function InputHelp({input, children}: {input: InputSpecification, children: any}) {
  return (
    <div className='input-help'>
      <div>{children}</div>
      {input.doc ? (
        <div className='help-icon-wrapper' title={input.doc}>
          <Icon icon='question-circle' className='help-icon'/>
        </div>
      ) : null}
    </div>
  );
}
