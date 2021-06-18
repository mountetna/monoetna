// Input component that takes nested object
//   and shows the keys one level at a time.
// Returns the last "Leaf" that the user selects.
import React, {useState, useEffect, useMemo, useCallback} from 'react';
import * as _ from 'lodash';

import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';
import {InputBackendComponent} from './input_types';
import TextInput from 'etna-js/components/inputs/text_input';

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

const ScatterPlotlyInput: InputBackendComponent = ({
  input,
  onChange
}) => {
  const {data, name} = input;

  const options: {[k: string]: string} = useMemo(() => {
    if (data) {
      return Object.values(data['data-options']);
    };

    return {};
  }, [data]);

  const defaults: {[k: string]: any} = {
    //'x_by': null,
    //'y_by': null,
    //'color_by': null,
    //'size': 5,
    //'color_order': 'increasing',
    //'order_when_continuous_color': false,
    'plot_title': 'make',
    'legend_title': 'make',
    'xlab': 'make',
    'ylab': 'make'
  }

  const updateValue = (newValue: string | number |  boolean, key: string, prevValues = input.value) => {
    prevValues[key] = newValue;
    onChange(name, prevValues);
  };

  const string_input = (key: string, label: string) => {
    <TextInput
      key={key}
      header={label}
      value={input.value[key]}
      onChange={(newValue: string) => updateValue(newValue, key)}
    />
  }

  /*
  const dropdown_input

  const checkbox_input

  const slider_input
  
  const component_use: (key: string, value: any) = {
    //'x_by': dropdown_input('X by', options),
    //'y_by': dropdown_input('Y by', options),
    //'color_by': dropdown_input('Color by', options),
    //'size': slider_input('Point Size', 0.1, 20),
    //'color_order': dropdown_input('Color Order', ['increasing', 'decreasing', 'unordered']),
    //'order_when_continuous_color': checkbox_input('Order point plotting when color is continuous?'),
    'plot_title': string_input('Plot Title', value),
    'legend_title': string_input('Legend Title', value),
    'xlab': string_input('X-Axis Title', value),
    'ylab': string_input('Y-Axis Title', value)
  }
  */
  const remove_hidden = (values: {[k: string]: any}, hide: [] | null) => {
    if (hide == null) {
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

  useEffect(() => {
    // Set all key's values as the initially given values, only if no input.value provided. (Initialization)
    // Otherwise set to mirror input.value
    if (null == input.value) {
      onChange(name, remove_hidden(defaults, data['hidden']));
    }
  }, [data, input.value, name, onChange]);

  if (null == input.value) return null;

  return (
    <div>
      {string_input('plot_title', 'Plot Title')}
      {string_input('legend_title', 'Legend Title')}
      {string_input('xlab', 'X-Axis Title')}
      {string_input('ylab', 'Y-Axis Title')}
    </div>
  );

};

export default ScatterPlotlyInput;
