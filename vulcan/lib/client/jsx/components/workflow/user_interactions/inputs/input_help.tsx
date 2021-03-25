import React from 'react';
import Icon from 'etna-js/components/icon';
import {InputSpecification} from "./input_types";

export default function InputHelp({input, children}: {input: InputSpecification, children: any}) {
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
