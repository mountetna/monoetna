import React, {useState, useCallback} from 'react';

import TextField from '@material-ui/core/TextField';

export function ShrinkingLabelTextField(props: any) {
  const [error, setError] = useState(false);

  const handleValidation = useCallback(
    (e: React.ChangeEvent<any>) => {
      if (!props.pattern) {
        props.onChange(e);
      } else {
        const input = e.target.value;

        if (!input.match(new RegExp(props.pattern))) {
          setError(true);
        } else {
          setError(false);
          props.onChange(e);
        }
      }
    },
    [props]
  );

  return (
    <TextField
      {...props}
      error={error}
      InputLabelProps={{
        shrink: props.value ? true : false
      }}
      inputProps={{
        pattern: props.pattern
      }}
      onChange={handleValidation}
    />
  );
}