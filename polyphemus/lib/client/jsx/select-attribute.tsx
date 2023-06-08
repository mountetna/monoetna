import React, { useContext } from 'react';
import Select from '@material-ui/core/Select';
import MenuItem from '@material-ui/core/MenuItem';
import {MagmaContext} from 'etna-js/contexts/magma-context';

const SelectAttribute = ({value, update, modelName, filter}:{
  value: string;
  update: Function;
  modelName: string;
  filter: Function;
}) => {
  const {models} = useContext(MagmaContext);

  const attribute_names = Object.values(
    models[modelName]?.template?.attributes || {}
  )
    .filter((a: any) => !a.hidden && (!filter || filter(a)))
    .map((a: any) => a.name)
    .sort();

  return <Select
    displayEmpty
    value={value}
    onChange={(e) => update(e.target.value as string)}
  >
    <MenuItem value='' disabled>
      <em>Select {modelName} attribute</em>
    </MenuItem>
    {attribute_names.map(att_name => (
      <MenuItem key={att_name} value={att_name}>
        {att_name}
      </MenuItem>
    ))}
  </Select>;
};

export default SelectAttribute;
