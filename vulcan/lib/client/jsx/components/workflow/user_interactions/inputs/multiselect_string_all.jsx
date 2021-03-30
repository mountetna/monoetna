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
      <MultiselectStringInput
        onAll={handleAllInputs}
        onClear={handleClearInputs}
        onChange={onChange}
        input={patchedInput} />
    </div>
  );
}
