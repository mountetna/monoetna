// Input component that takes nested object
//   and shows the keys one level at a time.
// Returns the last "Leaf" that the user selects.
import React, {useCallback, useMemo, useState} from 'react';
import * as _ from 'lodash';

import {DataEnvelope, nulled_vals, WithInputParams} from './input_types';
import { useSetsDefault } from './useSetsDefault';
import { Maybe, maybeOfNullable, some, withDefault } from '../../../../selectors/maybe';
import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';
import { Button, Slider } from '@material-ui/core';
import { pick } from 'lodash';
import NestedSelectAutocompleteInput from './nested_select_autocomplete';
import StringInput from './string';
import BooleanInput from './boolean';

/*
This input is closely tied to archimedes/functions/plotting/scatter_plotly.

Major design notes:
- It will have a set of advanced options that are shown/hidden via a toggle button
- It will also allow workflow-designers to hide controls for any inputs that are set internally by the workflow. (Ex: x_by & y_by for a UMAP)

Input structure:
  (Minimal = Either of 'data_frame' or 'data_options')
  'data_options': dictionary of options for '*_by' (ex: 'x_by', 'color_by') inputs where keys are the primary options and None vs [] values indicates what secondary options there may be if there are indeed secondary options (such as for genes!)
  'data_frame': A hash whose keys will be used as options if 'data_options' is not set
  'hide': An  array (python list). The names of scatter_plotly inputs which should not be shown in the current widget render.  Hidden inputs should be set elsewhere.  Perhaps they are hardset by the workflow so a user has no choice

Output Structure:
  dictionary of scatter_plotly inputs-name (key) + value pairs.
*/

export default function scatter_plotly({
  data, onChange, value
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  return VisualizationUI({data, onChange, value}, "scatterPlot")
}

function VisualizationUI({
  data, onChange, ...props
}: WithInputParams<{}, DataEnvelope<any>, any>, plotType: string) {
  const hide = useMemo(() => data && data['hide'], [data]);
  const defaultValue = useMemo(() => whichDefaults('scatter_plot', hide), [hide]);
  const value = useSetsDefault(defaultValue, props.value, onChange);
  const [showAdvanced, setShowAdvanced] = useState(false);

  const options: DataEnvelope<string> = useMemo(() => {
    if (data == null) return {};

    return data['data_options'] || nulled_vals(data['data_frame']);
  }, [data]);

  const updateValue = (newValue: any, key: string, prevValues = value) => {
    prevValues[key] = newValue;
    onChange(some(prevValues));
  };

  const extra_inputs = useExtraInputs(options);

  // Component set constructor
  const component_use = (key: string, value: any, extra_inputs: any) => {
    
    const comp_use: Function = comps[key]
    return(
      comp_use(key, updateValue, value, ...extra_inputs)
    )
  };

  // Advanced Options Button
  const base = useMemo(() => {
    return Object.keys(remove_hidden(
      {'x_by': 0, 'y_by': 0, 'color_by': 0}, hide
    ))
  }, [hide])
  
  const showHide: string = useMemo(() => {
    return (showAdvanced ? 'Hide' : 'Show')
  }, [showAdvanced])

  function toggleAdvanced() {
    setShowAdvanced(!showAdvanced)
  }

  const shownValues = useMemo(() => {
    return showAdvanced ? value : pick(value, ...base)
  }, [value, showAdvanced])
  
  return (
    <div>
      {Object.entries(shownValues).map(([key, val]) => {
        return component_use(key, val, extra_inputs[key])
      })}
      <div>
        <Button
          variant="contained"
          color="primary"
          onClick={() => {toggleAdvanced()}}
          >
          {showHide} Advanced Options
        </Button>
      </div>
    </div>
  );

};

const remove_hidden = (vals: DataEnvelope<any>, hide: string[] | null | undefined) => {
  
  let values = {...vals};

  if (hide == null || hide.length === 0) {
    return values;
  }
  
  const keys = Object.keys(values)
  for (let ind = 0; ind < keys.length; ind++) {
    if (hide.includes(keys[ind])) delete values[keys[ind]];
  }
  return values;
};

function val_wrap(v: any): DataEnvelope<typeof v> {
  return {'a': v}
}

function key_wrap(k: string[]) {
  let de: DataEnvelope<string> = {};
  for (let ind = 0; ind < k.length; ind++) {
    de[k[ind]]="0";
  }
  return de;
}

const input_sets: DataEnvelope<DataEnvelope<string[]>> = {
  scatter_plot: {
    main: ["x_by", "y_by", "color_by"],
    adv: ['size', 'plot_title', 'legend_title', 'xlab', 'ylab', 'color_order', 'order_when_continuous_color']
  }
}

const defaults: DataEnvelope<any> = {
  'x_by': null,
  'y_by': null,
  'color_by': null,
  'size': 5,
  'plot_title': 'make',
  'legend_title': 'make',
  'xlab': 'make',
  'ylab': 'make',
  'color_order': 'increasing',
  'order_when_continuous_color': false
};

function whichDefaults(plotType: string, hide: string[] | null | undefined) {
  const inputs = input_sets[plotType]['main'].concat(input_sets[plotType]['adv'])

  let output = defaults;
  const keys = Object.keys(defaults);
  for (let ind = 0; ind < defaults.length; ind++) {
    if (!inputs.includes(keys[ind])) delete output[keys[ind]];
  }

  return remove_hidden(defaults, hide);
}

function useExtraInputs(options: DataEnvelope<string>) {
  const extra_inputs: DataEnvelope<any[]> = useMemo(() => {
    return {
      // label, then for any extras
      'plot_title': ['Plot Title'],
      'legend_title': ['Legend Title'],
      'xlab': ['X-Axis Title'],
      'ylab': ['Y-Axis Title'],
      'x_by': ['X-Axis Data', options],
      'y_by': ['Y-Axis Data', options],
      'color_by': ['Color Points By', options],
      'color_order': ['Point render order', null, ['increasing', 'decreasing', 'unordered']],
      'order_when_continuous_color': ['Follow selected render ordering when color is continuous?'],
      'size': ['Point Size', 0.1, 50]
    }
  }, [options]);

  return extra_inputs;
}

// Component Setups
const stringInput = (
  key: string = "filler", changeFxn: Function, value: string = "filler",
  label: string = 'hello') => {
    return (
      <StringInput
        key={key}
        label={label}
        data={val_wrap(value)}
        value={maybeOfNullable(value)}
        onChange={(newValue) => changeFxn(withDefault(newValue,'make'), key)}
      />
    )
  };

const nestableDropdownInput = (
  key: string = "filler", changeFxn: Function, value: string | null,
  label: string, options: DataEnvelope<null>) => {
    
    return(
      <NestedSelectAutocompleteInput
        key={key}
        label={label}
        data={{options}} 
        value={maybeOfNullable(value)}
        onChange={(val) => changeFxn(withDefault(val, null), key)}
      />
    )
  }

const checkboxInput = (
  key: string = "filler", changeFxn: Function, value: boolean = false,
  label: string) => {

    return(
      <BooleanInput
        key={key}
        label={label}
        value={maybeOfNullable(value)}
        data={val_wrap(value)}
        onChange={() => changeFxn(!value, key)}
      />
    )
  }

const dropdownInput = (
  key: string = "filler", changeFxn: Function, value: string | null,
  label: string, defaultValue: string | null = null, options: string[]) => {
    
    return(
      <div key={key}>
        {label}
        <DropdownAutocomplete
          list={options}
          defaultValue={defaultValue}
          value={value}
          onSelect={(val: string) => changeFxn(val, key)}
        />
      </div>
    )
  }

const sliderInput = (
  key: string = "filler", changeFxn: Function, value: number,
  label: string, min: number = 0.1, max: number = 20) => {

    return(
        <div key={key}>
          {label}
          <Slider
            value={value}
            onChange={(event, newValue) => changeFxn(newValue as number, key)}
            min={min}
            max={max}
            valueLabelDisplay="auto"
          />
        </div>
    )
  }

const comps: DataEnvelope<Function> = {
  'plot_title': stringInput,
  'legend_title': stringInput,
  'xlab': stringInput,
  'ylab': stringInput,
  'x_by': nestableDropdownInput,
  'y_by': nestableDropdownInput,
  'color_by': nestableDropdownInput,
  'color_order': dropdownInput,
  'order_when_continuous_color': checkboxInput,
  'size': sliderInput
}