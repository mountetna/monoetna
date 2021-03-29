import React, {useEffect, useState} from 'react';

import MultiselectStringInput from './multiselect_string';

export default function MultiselectStringAllInput({input, onChange}) {
  const [patchedInput, setPatchedInput] = useState(null);

  if (!input || !onChange) return null;

  useEffect(() => {
    setPatchedInput({...input});
  }, [input]);

  function handleAllInputs() {
    let allOptions = input.options || [];
    setPatchedInput({
      ...input,
      default: allOptions
    });
    onChange(input.name, allOptions);
  }

  function handleClearInputs() {
    setPatchedInput({
      ...input,
      default: []
    });
    onChange(input.name, null);
  }

  return (
    <div>
      <div className='all-clear-header'>
        <div className='action' onClick={handleAllInputs}>
          All
        </div>
        <div className='action' onClick={handleClearInputs}>
          Clear
        </div>
      </div>
      <MultiselectStringInput onChange={onChange} input={patchedInput} />
    </div>
  );
}
