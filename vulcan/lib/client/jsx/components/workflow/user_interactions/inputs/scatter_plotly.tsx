// Input component that takes nested object
//   and shows the keys one level at a time.
// Returns the last "Leaf" that the user selects.
import React, {useState, useEffect, useMemo, useCallback} from 'react';
import * as _ from 'lodash';

import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';
import {InputBackendComponent} from './input_types';
import TextInput from 'etna-js/components/inputs/text_input';
import { JsxElement } from 'typescript';

/*
This input is closely tied to archimedes/functions/plotting/scatter_plotly.

Major design notes:
- It will have a set of advanced options that are shown/hidden via a toggle button
- It will also allow workflow-designers to hide controls for any inputs that are set internally by the workflow. (Ex: x_by & y_by for a UMAP)

Input structure:
  'data-options': dictionary of options for '*_by' (ex: 'x_by', 'color_by') inputs where keys are the primary options and None vs [] values indicates what secondary options there may be if there are indeed secondary options (such as for genes!)
  'hidden': list or null. The names of scatter_plotly inputs which should not be shown in the current widget render.  Perhaps these are hard-set by the workflow and a user has no choice.

Output Structure:
  dictionary of scatter_plotly inputs-name (key) + value pairs.
*/

const ScatterPlotly: InputBackendComponent = ({
  input,
  onChange
}) => {
  const {data, name} = input;

  const options: {[k: string]: string} = useMemo(() => {
    if (data == null) return {};
    
    let given = Object.values(data['data_options']);
    let cols = Object.keys(data['data_frame'])

    if (given == null && cols != null) {
      return cols
    }

    return given;
  }, [data]);

  const updateValue = (newValue: string | number | boolean, key: string, prevValues = input.value) => {
    prevValues[key] = newValue;
    console.log(input);
    onChange(name, prevValues);
  };

  const string_input = (key: string = "filler", value: string | number | boolean = "filler", label: string = 'hello') => {
    return (
      <TextInput
        key={key}
        header={label}
        value={value}
        onChange={(newValue: string) => updateValue(newValue, key)}
      />
    )
  };

  /*
  const dropdown_input

  const checkbox_input

  const slider_input
  */

  type fxn = (...args: any[]) => any;

  const defaults: {[k: string]: any} = {
    'plot_title': 'make',
    'legend_title': 'make',
    'xlab': 'make',
    'ylab': 'make'
  };
  
  const component_use = (key: string, value: any, extra_inputs: any) => {
    const comps: {[k: string]: fxn} = {
      'plot_title': string_input,
      'legend_title': string_input,
      'xlab': string_input,
      'ylab': string_input}
    
    const comp_use: fxn = comps[key]
    
    return(
      comp_use(key, value, ...extra_inputs)
    )
  };

  const extra_inputs: {[k: string]: any} = {
    // label, any extras
    'plot_title': ['Plot Title'],
    'legend_title': ['Legend Title'],
    'xlab': ['X-Axis Title'],
    'ylab': ['Y-Axis Title']
  };

  /*
  const remove_hidden = (values: {[k: string]: any}, hide: []) => {
    if (hide == null || hide === []) {
      return values;
    };
    
    let key;
    for (key in hide) {
      if (key in values.keys()) {
        delete values[key];
      };
    };
    return values;
  };
  */
  
      /*
    let hide;
    if (data != null) {
      hide = data['hidden'];
    } else {
      hide = [];
    }

    values = defaults;

    if (null == input.value) {
      onChange(name, remove_hidden(defaults, hide));
    }
    */
  useEffect(() => {
    if (null == input.value) {
      onChange(name, defaults);
    }
  }, [defaults, input.value, name, onChange]);

  if (null == input.value) return null;

  return (
    <div>
      {Object.entries(input.value).map(([key, val]) => {
        return (
          component_use(key, val, extra_inputs[key])
        );
      })}
    </div>
  );

};

export default ScatterPlotly;
