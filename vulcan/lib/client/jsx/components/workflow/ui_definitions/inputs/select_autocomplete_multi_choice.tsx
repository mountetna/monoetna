import React, {useCallback, useEffect, useMemo, useState} from 'react';
import {WithInputParams} from '../input_types';
import {maybeOfNullable, some, withDefault} from '../../../../selectors/maybe';
import {flattenStringOptions, StringOptions} from '../monoids';
import {useMemoized} from '../../../../selectors/workflow_selectors';
import {useSetsDefault} from '../useSetsDefault';
import TextField from '@material-ui/core/TextField';
import Autocomplete from '@material-ui/lab/Autocomplete';
import {useAsyncCallback} from 'etna-js/utils/cancellable_helpers';
import { pullRecommendation } from './components/select_autocomplete'

// Compared to select_autocomplete, also filter selected options
function filterOptions(query: string, opts: string[], value: string[]) {
  return opts.filter((o) => {
    return query == null ? true : !value.includes(o) && o.toLowerCase().indexOf(query.toLowerCase()) > -1;
  });
}
const getOptionsAsync = (query: string, opts: string[], value: string[]) => {
  return new Promise((resolve) => {
    setTimeout(() => {
      resolve(filterOptions(query, opts, value));
    }, 3000);
  });
};

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

export default function SelectAutocompleteMultiPickInput({
  data,
  label,
  minWidth,
  maxOptions = 100,
  disableClearable = true,
  disabled = false,
  placeholder,
  onChangeOverride,
  testId,
  onChange,
  ...props
}: WithInputParams<
  {
    label?: string;
    minWidth?: number;
    disableClearable?: boolean;
    disabled?: boolean;
    placeholder?: string;
    maxOptions?: number;
    onChangeOverride?: (event: any, e: string[]) => void;
    testId?: string
  },
  string[],
  StringOptions
>) {
  /*
  Creates a searchable dropdown selection input box from concatenated values of the 'data' hash.
  Special Case: If any data key is "recommendation", a line of text will display the values of this recommendation to the user.
  */
  const [data_use, suggestion] = useMemoized(pullRecommendation, data);

  const options_in = useMemoized(flattenStringOptions, data_use);
  const value = useSetsDefault([], props.value, onChange) as string[];
  const disp_label = useMemo(() => {
    return suggestion ? suggestion : label;
  }, [suggestion, label]);
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
      getOptionsDelayed(query, options_in, value, (filteredOptions: string[]) => {
        setLoadingOptions(false);
        setOptions(filteredOptions);
      });
    } else {
      // console.log('calculating options - fast')
      setOptions(filterOptions(query, options_in, value));
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

  const onChangeAction = useCallback(
    (event: any, e: string[]) => {
      onChangeOverride
        ? onChangeOverride(event, e)
        : onChange(maybeOfNullable(e));
    },
    [onChangeOverride, onChange]
  );

  // ToDo: Replace with stubbing options_in + disabled + placeholder text to explain
  if (0 === Object.keys(data_use || {}).length) return null;

  return (
    <Autocomplete
      key={label}
      disableClearable={disableClearable}
      disabled={disabled}
      multiple
      clearOnBlur={true}
      options={options_in}
      filterOptions={(x: string[]) => options.display}
      loading={loadingOptions}
      value={value==null ? [] : value}
      onChange={onChangeAction}
      inputValue={inputState}
      onInputChange={(event: any, newInputState: string) => {
        setInputState(newInputState);
      }}
      id={testId ? `${testId}-full` : 'multi-pick-input-full'}
      style={{minWidth: minWidth, paddingTop: disp_label ? 8 : 0}}
      renderInput={(params: any) => (
        <TextField
          {...params}
          helperText={helperText}
          error={value==null ? true : value.includes(inputState)}
          label={disp_label}
          placeholder={placeholder}
          fullWidth
          size='small'
          InputLabelProps={{shrink: true}}
          inputProps={{
            'data-testid': testId ? testId : 'multi-pick-input',
            'id': testId ? testId : 'multi-pick-input',
          }}
        />
      )}
    />
  );
}
