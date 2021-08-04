// Input component that takes nested object
//   and shows the keys one level at a time.
// Returns the last "Leaf" that the user selects.
import React, {useCallback, useMemo} from 'react';
import * as _ from 'lodash';

import {DataEnvelope, get_keys, get_vals, nulled_vals, WithInputParams} from './input_types';
import TextInput from 'etna-js/components/inputs/text_input';
import { useSetsDefault } from './useSetsDefault';
import { applySome, some } from '../../../../selectors/maybe';
import DropdownAutocomplete from 'etna-js/components/inputs/dropdown_autocomplete';
import { useEffect } from 'react';
import MultiselectStringInput from './multiselect_string';
import { Button, Slider } from '@material-ui/core';
import { pick } from 'lodash';

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

const remove_hidden = (values: DataEnvelope<any>, hide: string[] | null | undefined) => {
  if (hide == null || hide.length === 0) {
    return values;
  }
  
  const keys = Object.keys(values)
  for (let ind = 0; ind < keys.length; ind++) {
    if (hide.includes(keys[ind])) delete values[keys[ind]];
  }
  return values;
};

let showAdvanced = false;

export default function ScatterPlotly({
  data, onChange, ...props
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  const hide = useMemo(() => data && data['hide'], [data]);
  const defaultValue = useMemo(() => remove_hidden(defaults, hide), [hide]);
  const value = useSetsDefault(defaultValue, props.value, onChange);

  const options: DataEnvelope<string> = useMemo(() => {
    if (data == null) return {};

    return data['data_options'] || nulled_vals(data['data_frame']);
  }, [data]);

  const updateValue = (newValue: string | number | boolean, key: string, prevValues = value) => {
    prevValues[key] = newValue;
    onChange(some(prevValues));
  };

  // Component Setups
  const string_input = (
    key: string = "filler", value: string | number | boolean = "filler",
    label: string = 'hello') => {
      return (
        <TextInput
          key={key}
          header={label}
          value={value}
          onChange={(newValue: string) => updateValue(newValue, key)}
        />
      )
    };

  const dropdown_input = (
    key: string = "filler", value: string | number,
    label: string, options: DataEnvelope<null>) => {
      
      //console.log(key, label);
      //console.log(key, options)

      // Need to make this work without the `get_keys` similar to nested_select_autocomplete.tsx component
      return(
        <div key={key}>
          {label}
          <DropdownAutocomplete
            onSelect={(newValue: string) => updateValue(newValue, key)}
            list={get_keys(options)}
            value={value}
          />
        </div>
        
      )

      /*
      const handleSelect = useCallback(
        (value: string | null, depth: number) => {
          // User has not selected something...perhaps
          //   still typing?
          if (null == value) return;
          if (options == null) return;
    
          // figure out where the next selection comes from, check if we're picking a leaf node.
          const updatedOptions = options.slice(0, depth);
          updatedOptions.push(value);
          setPath(updatedPath);
    
          if (getOptions(updatedPath, allOptions) == null) {
            // If we are updating a leaf
            onChange(some(value));
          } else {
            // Otherwise a leaf has not been selected.
            onChange(null);
          }
        },
        [allOptions, path, onChange]
      );
    
      return (
        <div>
          <div>
            {
              <DropdownAutocomplete
                key={key}
                onSelect={(e: string | null) => {
                  handleSelect(e, index);
                }}
                list={options}
                value={value}
              />
            }
          </div>
          <div>
            <LeafOptions
              options={getOptions(path, allOptions)}
              depth={path.length}
              handleSelect={handleSelect}
            />
          </div>
        </div>
      );
      */
    }

  const checkbox_input = (
    key: string = "filler", value: boolean = false,
    label: string) => {

      return(
          <div key={key}>
            <text>{label}</text>
            <input id={key} type='checkbox'
                  checked={value}
                  onChange={() => updateValue(!value, key)} />
          </div>
      )
    }

  const slider_input = (
    key: string = "filler", value: number,
    label: string, min: number = 0.1, max: number = 20) => {

      return(
          <div key={key}>
            <text>{label}</text>
            <Slider
              value={value}
              onChange={(event, newValue) => updateValue(newValue as number, key)}
              min={min}
              max={max}
              valueLabelDisplay="auto"
            />
          </div>
      )
    }

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
      'color_order': ['Point render order', {'increasing': null, 'decreasing': null, 'unordered': null}],
      'order_when_continuous_color': ['Follow selected render ordering when color is continuous?'],
      'size': ['Point Size', 0.1, 50]
    }
  }, [options]);
  //console.log("extra_inputs:", extra_inputs)

  // Component set constructor
  const component_use = (key: string, value: any, extra_inputs: any) => {
    const comps: DataEnvelope<Function> = {
      'plot_title': string_input,
      'legend_title': string_input,
      'xlab': string_input,
      'ylab': string_input,
      'x_by': dropdown_input,
      'y_by': dropdown_input,
      'color_by': dropdown_input,
      'color_order': dropdown_input,
      'order_when_continuous_color': checkbox_input,
      'size': slider_input
    }
    
    const comp_use: Function = comps[key]
    //console.log("component:", key, value, ...extra_inputs)
    return(
      comp_use(key, value, ...extra_inputs)
    )
  };

  // Advanced Options Button
  const base = useMemo(() => {
    return get_keys(remove_hidden(
      {'x_by': 0, 'y_by': 0, 'color_by': 0}, hide
    ))
  }, [hide])
  
  const showHide: string = useMemo(() => {
    return (showAdvanced ? 'Hide' : 'Show')
  }, [showAdvanced])

  function toggleAdvanced() {
    showAdvanced = !showAdvanced
    onChange(some(value))
  }

  const shownValues = useMemo(() => {
    return showAdvanced ? value : pick(value, ...base)
  }, [value, showAdvanced])

  //console.log('values:', value)
  //console.log('true values:', props.value)
  //console.log('shownValues', shownValues)
  
  return (
    <div>
      {Object.entries(shownValues).map(([key, val]) => {
        return component_use(key, val, extra_inputs[key])
      })}
      <Button
        variant="contained"
        color="primary"
        onClick={() => {toggleAdvanced()}}
        >
         {showHide} Advanced Options
      </Button>
    </div>
  );

};

