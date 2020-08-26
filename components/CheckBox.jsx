import React from 'react';
require('./CheckBox.css');

export function CheckBox(props) {
  return <input  type='checkbox' className='etna-checkbox' {...props}/>
}