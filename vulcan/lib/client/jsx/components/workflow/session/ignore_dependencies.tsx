import React from 'react';
import {Checkbox, FormControlLabel} from '@material-ui/core';

export default function IgnoreDependencies({
  onChange,
  value
}: {
  onChange: () => void;
  value: boolean;
}) {
  const inner = (
    <Checkbox
      onChange={onChange}
      checked={value}
      inputProps={{'aria-label': 'controlled'}}
    />
  );

  return <FormControlLabel control={inner} label='Ignore dependencies' />;
}
