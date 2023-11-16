import React, {useEffect, useState} from 'react';
import TextField from '@material-ui/core/TextField';
import Autocomplete from '@material-ui/lab/Autocomplete';
import {useAsyncCallback} from 'etna-js/utils/cancellable_helpers';

export function determineHelperText(
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

function filterOptions(query: string, opts: string[], values: string[]) {
  return opts.filter((o) => {
    return query == null ? !values.includes(o) : !values.includes(o) && o.indexOf(query) > -1;
  });
}

const getOptionsAsync = (query: string, opts: string[], values: string[]) => {
  return new Promise((resolve) => {
    setTimeout(() => {
      resolve(filterOptions(query, opts, values));
    }, 3000);
  });
};

export default function DropdownAutocompleteMultiPickInput({
  optionSet,
  value,
  label,
  placeholder,
  minWidth,
  maxOptions = 100,
  disableClearable = true,
  disabled = false,
  enforceError = false,
  onChange,
  ...props
}: {
  optionSet: string[];
  value: string[]
  label?: string;
  placeholder?: string;
  minWidth?: number;
  disableClearable?: boolean;
  disabled?: boolean;
  enforceError?: boolean;
  maxOptions?: number;
  onChange?: (event: any, e: string[]) => void;
}) {
  /*
  Creates a searchable dropdown selection input box from optionSet values.
  The options display is a scrollable box that is capped at the first 100 options by default to keep the rendering light-weight.
  */

  const [loadingOptions, setLoadingOptions] = useState(false);
  const [helperText, setHelperText] = useState(undefined as string | undefined);
  const [inputState, setInputState] = useState('');
  const [options, setOptionsFull] = useState({
    filtered: optionSet,
    display: [...optionSet].splice(0, maxOptions)
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
    const options = yield getOptionsAsync(text, [...options_in], value);
    callback(options);
  },
  []);

  useEffect(() => {
    // Shown from all options when user has made their selection (current text = current value)
    const query = value==null || value.includes(inputState) ? '' : inputState;

    if (optionSet.length > 1000) {
      setLoadingOptions(true);
      // console.log('calculating options - slow')
      getOptionsDelayed(query, optionSet, (filteredOptions: string[]) => {
        setLoadingOptions(false);
        setOptions(filteredOptions);
      });
    } else {
      // console.log('calculating options - fast')
      setOptions(filterOptions(query, optionSet, value));
    }
  }, [inputState, getOptionsDelayed, optionSet, value]);

  useEffect(() => {
    determineHelperText(
      setHelperText,
      options.filtered,
      optionSet,
      loadingOptions,
      maxOptions,
      inputState,
      value
    );
  }, [
    options,
    optionSet,
    inputState,
    value,
    maxOptions,
    loadingOptions,
    setHelperText
  ]);

  return (
    <Autocomplete
      key={label}
      disableClearable={disableClearable}
      disabled={disabled}
      multiple
      clearOnBlur={true}
      options={optionSet}
      filterOptions={(x: string[]) => options.display}
      loading={loadingOptions}
      value={value}
      onChange={onChange}
      inputValue={inputState}
      onInputChange={(event: any, newInputState: string) => {
        setInputState(newInputState);
      }}
      style={{minWidth: minWidth, paddingTop: label ? 8 : 0}}
      renderInput={(params: any) => (
        <TextField
          {...params}
          helperText={helperText}
          error={enforceError || value.includes(inputState)}
          label={label}
          placeholder={value.length > 0 ? undefined : placeholder}
          size='small'
          InputLabelProps={{shrink: true}}
        />
      )}
    />
  );
}