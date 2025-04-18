import React, {useCallback, useEffect, useMemo, useState} from 'react';
import TextField from '@material-ui/core/TextField';
import Autocomplete from '@material-ui/lab/Autocomplete';
import {useAsyncCallback} from 'etna-js/utils/cancellable_helpers';
import { PieceBaseInputs } from './user_input_pieces';

// Compared to select_autocomplete, also filter selected options
function valueUse(value: string[] | null) {
  return value==null ? [] : value;
};
function filterOptions(query: string, opts: string[], value: string[]) {
  return opts.filter((o) => {
    return query == null ? true : !value.includes(o) && o.toLowerCase().indexOf(query.toLowerCase()) > -1;
  });
};
const getOptionsAsync = (query: string, opts: string[], value: string[]) => {
  return new Promise((resolve) => {
    setTimeout(() => {
      resolve(filterOptions(query, opts, value));
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
    saved: string[] | null
  ) {
    if (loadingOptions) {
      return;
    } else if (saved == null || saved.includes(typed)) {
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

export interface DropdownMultiChoicePieceInputs extends PieceBaseInputs<string[] | null> {
  options_in: string[],
  minWidth?: number,
  maxOptions?: number,
  disableClearable?: boolean,
  disabled?: boolean,
  placeholder?: string,
  testId?: string
}

export function DropdownMultiChoicePiece(
  key: DropdownMultiChoicePieceInputs['name'],
  changeFxn: DropdownMultiChoicePieceInputs['changeFxn'],
  value: DropdownMultiChoicePieceInputs['value'],
  label: DropdownMultiChoicePieceInputs['label'] = '',
  options: DropdownMultiChoicePieceInputs['options_in'],
  minWidth: DropdownMultiChoicePieceInputs['minWidth'] = 200,
  disabled: DropdownMultiChoicePieceInputs['disabled'] = false,
): React.ReactElement | null {
  return <DropdownMultiChoicePieceRct
    name={key}
    changeFxn={changeFxn}
    value={value}
    label={label}
    options_in={options}
    minWidth={minWidth}
    disabled={disabled}
  />
}

export function DropdownMultiChoicePieceRct({
  name,
  changeFxn,
  value,
  label = '',
  options_in,
  minWidth,
  maxOptions = 200,
  disableClearable = true,
  disabled = false,
  placeholder = undefined,
  testId = undefined
}: DropdownMultiChoicePieceInputs): React.ReactElement | null {

  const [loadingOptions, setLoadingOptions] = useState(false);
  const [helperText, setHelperText] = useState(undefined as string | undefined);
  const [inputState, setInputState] = useState('');
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
    value: string[],
    callback: Function
  ) {
    const options = yield getOptionsAsync(text, [...options_in], value);
    callback(options);
  },
  []);

  useEffect(() => {
    // Shown from all options when user has made their selection (current text = current value)
    const query = value==null || value.includes(inputState) ? '' : inputState;

    if (options_in.length > 1000) {
      setLoadingOptions(true);
      // console.log('calculating options - slow')
      getOptionsDelayed(query, options_in, valueUse(value), (filteredOptions: string[]) => {
        setLoadingOptions(false);
        setOptions(filteredOptions);
      });
    } else {
      // console.log('calculating options - fast')
      setOptions(filterOptions(query, options_in, valueUse(value)));
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

  return (
    <Autocomplete
      key={name}
      disableClearable={disableClearable}
      disabled={disabled}
      multiple
      clearOnBlur={true}
      options={options_in}
      filterOptions={(x: string[]) => options.display}
      loading={loadingOptions}
      value={valueUse(value)}
      onChange={(event: any, e: string[]) => changeFxn(e, name)}
      inputValue={inputState}
      onInputChange={(event: any, newInputState: string) => {
        setInputState(newInputState);
      }}
      id={!!testId ? `${testId}` : name}
      style={{minWidth: minWidth}}
      renderInput={(params: any) => (
        <TextField
          {...params}
          helperText={helperText}
          error={value==null ? true : value.includes(inputState)}
          label={label}
          placeholder={placeholder}
          fullWidth
          size='small'
          InputLabelProps={{shrink: true}}
          id={!!testId ? `${testId}` : undefined}
          data-testid={!!testId ? `${testId}` : undefined}
        />
      )}
    />
  );
}
