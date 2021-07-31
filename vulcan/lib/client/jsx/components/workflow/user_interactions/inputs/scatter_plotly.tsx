// Input component that takes nested object
//   and shows the keys one level at a time.
// Returns the last "Leaf" that the user selects.
import React, {useMemo} from 'react';
import * as _ from 'lodash';

import {DataEnvelope, WithInputParams} from './input_types';
import TextInput from 'etna-js/components/inputs/text_input';
import { useSetsDefault } from './useSetsDefault';
import { some } from '../../../../selectors/maybe';
import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';

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

const defaults: {[k: string]: any} = {
  'x_by': null,
  'y_by': null,
  'color_by': null,
  'plot_title': 'make',
  'legend_title': 'make',
  'xlab': 'make',
  'ylab': 'make'
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

export default function ScatterPlotly({
  data, onChange, ...props
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  const value = useSetsDefault(defaults, props.value, onChange);

  const options: DataEnvelope<string> = useMemo(() => {
    if (data == null) return {};
    
    let given = Object.values(data['data_options']);
    let cols = Object.keys(data['data_frame'])

    return cols || given;
  }, [data]);

  const updateValue = (newValue: string | number | boolean, key: string, prevValues = value) => {
    prevValues[key] = newValue;
    onChange(some(prevValues));
  };

  // Component Setups
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

  const dropdown_input = (key: string = "filler", value: string | number, label: string) => {
    return(
      <div>
        <text>{label}</text>
        <DropdownAutocomplete
          key={key}
          list={options}
          onSelect={(newValue: string) => updateValue(newValue, key)}
        />
      </div>
    )
  }
  /*
  const checkbox_input

  const slider_input
  */

  const extra_inputs: {[k: string]: any} = {
    // label, any extras
    'plot_title': ['Plot Title'],
    'legend_title': ['Legend Title'],
    'xlab': ['X-Axis Title'],
    'ylab': ['Y-Axis Title'],
    'x_by': ['X-Axis Data'],
    'y_by': ['Y-Axis Data'],
    'color_by': ['Point Color Data']
  };

  // Component set constructor
  const component_use = (key: string, value: any, extra_inputs: any) => {
    const comps: DataEnvelope<Function> = {
      'plot_title': string_input,
      'legend_title': string_input,
      'xlab': string_input,
      'ylab': string_input,
      'x_by': dropdown_input,
      'y_by': dropdown_input,
      'color_by': dropdown_input
    }
    
    const comp_use: Function = comps[key]
    
    return(
      comp_use(key, value, ...extra_inputs)
    )
  };

  return (
    <div>
      {Object.entries(value).map(([key, val]) => {
        return (
          component_use(key, val, extra_inputs[key])
        );
      })}
    </div>
  );

};
