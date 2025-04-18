import ListInput from 'etna-js/components/inputs/list_input';
import DropdownAutocompleteInput from 'etna-js/components/inputs/dropdown_autocomplete_wrapper';
import React from 'react';
import InputLabel from '@material-ui/core/InputLabel';
import { DataEnvelope } from '../../input_types';
import TextField from '@material-ui/core/TextField';
import { arrayLevels, DisabledTextbox, PieceBaseInputs } from './user_input_pieces';

interface MultiselectStringPieceInputs extends PieceBaseInputs<string[] | null> {
  options: string[];
  onAll?: string[];
  onClear?: string[];
  maxItems?: number;
}

export function MultiselectStringPieceRct({
  name,
  changeFxn,
  value,
  label,
  options,
  onAll = options,
  onClear = [],
  maxItems = 25
}: MultiselectStringPieceInputs): React.ReactElement {

  if (value != null && value.length > 0) {
    if (value.filter((val) => !(val==='' || options.includes(val))).length > 0) changeFxn(null, name);
  }

  return(
    <div key={name}>
      <InputLabel htmlFor={'multiselect-'+name} shrink>{label}</InputLabel>
      <ListInput
        key={'multiselect-'+name}
        placeholder='Select items from the list'
        className='link_text'
        values={value!=null ? value : []}
        itemInput={DropdownAutocompleteInput}
        list={options}
        onChange={(val: string[]) => changeFxn(val, name)}
        onAll={() => changeFxn(onAll, name)}
        onClear={() => changeFxn(onClear, name)}
        maxItems={maxItems}
      />
    </div>
  );
}

export function MultiselectStringPiece(
  key: MultiselectStringPieceInputs['name'],
  changeFxn: MultiselectStringPieceInputs['changeFxn'],
  value: MultiselectStringPieceInputs['value'],
  label: MultiselectStringPieceInputs['label'],
  options: MultiselectStringPieceInputs['options'],
  onAll?: MultiselectStringPieceInputs['onAll'],
  onClear?: MultiselectStringPieceInputs['onClear'],
  maxItems?: MultiselectStringPieceInputs['maxItems']
) {
  return <MultiselectStringPieceRct
    name={key}
    changeFxn={changeFxn}
    value={value}
    label={label}
    options={options}
    onAll={onAll}
    onClear={onClear}
    maxItems={maxItems}
  />
}

interface MultiselectStringAfterDataChoicePieceInputs extends PieceBaseInputs<string[] | null> {
  full_data: DataEnvelope<any[]>,
  data_target: string | null,
  data_target_label: string,
  discrete_data: string[]
}

export function MultiselectStringAfterDataChoicePieceRct({
  name,
  changeFxn,
  value,
  label,
  full_data,
  data_target,
  data_target_label,
  discrete_data
}: MultiselectStringAfterDataChoicePieceInputs): React.ReactElement {
  const canReorder = data_target != null && discrete_data.includes(data_target) && full_data != null
  const levels = canReorder ? arrayLevels(Object.values(full_data[data_target as string])) : null
  return(
    levels != null ? <MultiselectStringPieceRct
      name={name}
      changeFxn={changeFxn}
      value={value}
      label={label}
      options={levels}
    /> : <DisabledTextbox
      name={`multiselect-${name}`}
      label={label}
      text={`Awaiting ${data_target_label} choice`}
    />
  )
}

export function MultiselectStringAfterDataChoicePiece(
  key: MultiselectStringAfterDataChoicePieceInputs['name'],
  changeFxn: MultiselectStringAfterDataChoicePieceInputs['changeFxn'],
  value: MultiselectStringAfterDataChoicePieceInputs['value'],
  label: MultiselectStringAfterDataChoicePieceInputs['label'],
  full_data: MultiselectStringAfterDataChoicePieceInputs['full_data'],
  data_target: MultiselectStringAfterDataChoicePieceInputs['data_target'], // The name of a column/key of full_data which the user has chosen as the target of this ui piece, or null if not chosen yet.
  data_target_label: MultiselectStringAfterDataChoicePieceInputs['data_target_label'], // The label of the ui-piece where the user selects data_target 
  discrete_data: MultiselectStringAfterDataChoicePieceInputs['discrete_data']
) {
  return <MultiselectStringAfterDataChoicePieceRct
    name={key}
    changeFxn={changeFxn}
    value={value}
    label={label}
    full_data={full_data}
    data_target={data_target}
    data_target_label={data_target_label}
    discrete_data={discrete_data}
  />
}