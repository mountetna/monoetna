import React, {useCallback, useMemo, useState} from 'react';

import MultiselectStringInput, {getAllOptions} from './multiselect_string';
import {InputBackendComponent} from "./input_types";

const MultiselectStringAllInput: InputBackendComponent = ({input, onChange}) => {
  const options = useMemo(() => getAllOptions(input.data), [input.data]);
  const [defaults, setDefaults] = useState(input.default);
  const patchedInput = useMemo(() => ({ ...input, default: defaults }), [input, defaults]);

  const handleAllInputs = useCallback(() => {
    setDefaults(options);
    onChange(input.name, options);
  }, [input.name, options]);


  const handleClearInputs = useCallback(() => {
    setDefaults([]);
    onChange(input.name, null);
  }, [input.name]);

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


export default MultiselectStringAllInput;