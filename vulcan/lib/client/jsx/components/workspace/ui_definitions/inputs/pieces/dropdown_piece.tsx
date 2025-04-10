import React, {useEffect, useState} from 'react';
import TextField from '@material-ui/core/TextField';
import Autocomplete from '@material-ui/lab/Autocomplete';
import {useAsyncCallback} from 'etna-js/utils/cancellable_helpers';
import { PieceBaseInputs } from './user_input_pieces';

function dispValue(value: string | null) {
  return value == null ? '' : value;
}

export function filterOptions(query: string, opts: string[]) {
  return opts.filter((o) => {
    return query == null ? true : o.toLowerCase().indexOf(query.toLowerCase()) > -1;
  });
}

export const getOptionsAsync = (query: string, opts: string[]) => {
  return new Promise((resolve) => {
    setTimeout(() => {
      resolve(filterOptions(query, opts));
    }, 3000);
  });
};

function determineHelperText(
  setHelperText: React.Dispatch<React.SetStateAction<string | undefined>>,
  filteredOptions: string[],
  allOptions: string[],
  loadingOptions: boolean,
  maxOptions: number,
  typed: string,
  saved: string | null
) {
  if (loadingOptions) {
    return;
  } else if (saved == null || typed != saved) {
    if (filteredOptions.length > maxOptions) {
      setHelperText(
        allOptions.length + ' options, only ' + maxOptions + ' shown'
      );
    } else if (filteredOptions.length == 0) {
      setHelperText('no matching options');
    } else {
      setHelperText(undefined);
    }
  } else {
    setHelperText(undefined);
  }
}

interface DropdownPieceInputs extends PieceBaseInputs<string | null>{
  options_in: string[],
  minWidth?: number,
  maxOptions?: number,
  disableClearable?: boolean,
  disabled?: boolean
}

export default function DropdownPiece(
  key: DropdownPieceInputs['name'],
  changeFxn: DropdownPieceInputs['changeFxn'],
  value: DropdownPieceInputs['value'],
  label: DropdownPieceInputs['label'] = '',
  options: string[] | {[k: string]: string} | null,
  sorted: boolean = true,
  minWidth: DropdownPieceInputs['minWidth'] = 200,
  disabled: DropdownPieceInputs['disabled'] = false,
): React.ReactElement | null {
  if (options==null) return null;
  if (Array.isArray(options)) {
    return <DropdownPieceRct
      name={key}
      changeFxn={changeFxn}
      value={value}
      label={label}
      options_in={options}
      minWidth={minWidth}
      disabled={disabled}
    />
  }
  // options_in = Object with keys = label values to show the user; values = true option values that the output / saved state should hold.
  const labels = Object.keys(options)
  const toValue = (label: string | null | undefined) => label == null? null : options[label]
  const toLabel = (value: string | null) => {
    return (value == null) ? null : (Object.keys(options).find(key => options[key] === value)) as string
  };
  return(
    <DropdownPieceRct
      name={key}
      changeFxn={(v,k) => {changeFxn(toValue(v), k)}}
      label={label}
      value={toLabel(value)}
      options_in={labels}
      minWidth={minWidth}
    />
  );
}

export function DropdownPieceRct({
  name,
  changeFxn,
  value,
  label = '',
  options_in,
  minWidth,
  maxOptions = 200,
  disableClearable = true,
  disabled = false,
}: DropdownPieceInputs): React.ReactElement | null {

  if (!options_in) {
    console.log(`No Options Given to dropdown-piece named '${name}'`)
    return null
  }

  const [loadingOptions, setLoadingOptions] = useState(false);
  const [helperText, setHelperText] = useState(undefined as string | undefined);
  const [inputState, setInputState] = useState(dispValue(value));
  const [options, setOptionsFull] = useState({
    filtered: options_in,
    display: [...options_in].splice(0, maxOptions)
  });
  function setOptions(filteredOptions: string[]) {
    setOptionsFull({
      filtered: filteredOptions,
      display: [...filteredOptions].splice(0, maxOptions)
    });
  }
  const [getOptionsDelayed] = useAsyncCallback(function* (
    text: string,
    options_in: string[],
    callback: Function
  ) {
    const options = yield getOptionsAsync(text, [...options_in]);
    callback(options);
  },
  []);

  useEffect(() => {
    // Shown from all options when user has made their selection (current text = current value)
    const query = inputState != value ? inputState : '';

    if (options_in.length > 1000) {
      setLoadingOptions(true);
      // console.log('calculating options - slow')
      getOptionsDelayed(query, options_in, (filteredOptions: string[]) => {
        setLoadingOptions(false);
        setOptions(filteredOptions);
      });
    } else {
      // console.log('calculating options - fast')
      setOptions(filterOptions(query, options_in));
    }
  }, [inputState, getOptionsDelayed, options_in, value]);

  useEffect(() => {
    determineHelperText(
      setHelperText,
      options.filtered,
      options_in,
      loadingOptions,
      maxOptions,
      inputState,
      value
    );
  }, [
    options,
    options_in,
    inputState,
    value,
    maxOptions,
    loadingOptions,
    setHelperText
  ]);

  // ToDo: Replace with stubbing options_in + disabled + placeholder text to explain
  if (0 === options_in.length) return null;

  return (
    <div key={!!name ? name : 'dropdown'}>
      <Autocomplete
        disableClearable={disableClearable}
        disabled={disabled}
        clearOnBlur={true}
        options={options_in}
        filterOptions={(x: string[]) => options.display}
        loading={loadingOptions}
        value={value}
        onChange={(event: any, e: string | null) => {changeFxn(e, name)}}
        inputValue={inputState}
        onInputChange={(event: any, newInputState: string) => {setInputState(newInputState)}}
        style={{minWidth: minWidth, paddingTop: label!='' ? 8 : 0}}
        renderInput={(params: any) => (
          <TextField
            {...params}
            helperText={helperText}
            error={inputState != dispValue(value)}
            label={label}
            size='small'
            InputLabelProps={{shrink: true}}
          />
        )}
      />
    </div>
  );
}
