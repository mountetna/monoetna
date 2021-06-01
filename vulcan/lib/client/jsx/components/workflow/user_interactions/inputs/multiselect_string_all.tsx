import React, {useCallback, useMemo, useState} from 'react';

import MultiselectStringInput, {getAllOptions} from './multiselect_string';
import {InputBackendComponent} from "./input_types";

const MultiselectStringAllInput: InputBackendComponent = ({input, onChange}) => {
  const options = useMemo(() => getAllOptions(input.data), [input.data]);

  const handleAllInputs = useCallback(() => {
    onChange(input.name, options);
  }, [input.name, onChange, options]);


  const handleClearInputs = useCallback(() => {
    onChange(input.name, null);
  }, [input.name, onChange]);

  return (
    <div>
      <MultiselectStringInput
        onAll={handleAllInputs}
        onClear={handleClearInputs}
        onChange={onChange}
        input={input} />
    </div>
  );
}


export default MultiselectStringAllInput;