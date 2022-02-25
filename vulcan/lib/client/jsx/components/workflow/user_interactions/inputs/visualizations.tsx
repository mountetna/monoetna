// Input component that takes nested object
//   and shows the keys one level at a time.
// Returns the last "Leaf" that the user selects.
import React, {useMemo, useState} from 'react';
import * as _ from 'lodash';

import {DataEnvelope, WithInputParams} from './input_types';
import { useSetsDefault } from './useSetsDefault';
import { some } from '../../../../selectors/maybe';
import { Button, Slider } from '@material-ui/core';
import { pick } from 'lodash';
import { key_wrap, stringPiece, dropdownPiece, MultiselectPiece, checkboxPiece, sliderPiece } from './user_input_pieces';
import { subsetDataFramePiece } from './subsetDataFrame_piece';

/*
This UI is closely tied to archimedes/functions/plotting/*_plotly functions.

Major design notes:
- Organized around a set of pairs of input-names and associated component-setups.
- An 'input_sets' object then defines the set of inputs used for any given visualization type & a new exported function should be created for every such plot_type.
- Widgets have a set of advanced options that are shown/hidden via a toggle button. These inputs are listed as 'adv' inputs in the 'input_sets' definition.
- Widgets also allow workflow-designers to hide controls for any inputs that are set internally by the workflow. (Ex: x_by & y_by for a UMAP)

Input structure:
  'data_frame': A hash representing the data_frame that will be used to make a plot. keys = column names, values = data points. 
  (optional) 'preset': A hash where keys = input_names that should be hidden from the user & values = the preset value that should be used for that input. E.g. The x_by and y_by inputs are hardset within the umap workflow as '0' and '1', so a user has no choice here and should not be able to adjust those fields!

Output Structure:
  A hash (dict once in python) of input-name (key) + value pairs that can be splatted, along with the accompanying DataFrame, into a visualization funciton in archimedes.
*/

export function ScatterPlotly({
  data, onChange, value
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  return VisualizationUI({data, onChange, value}, "scatter_plot")
}

export function BarPlotly({
  data, onChange, value
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  return VisualizationUI({data, onChange, value}, "bar_plot")
}

export function YPlotly({
  data, onChange, value
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  return VisualizationUI({data, onChange, value}, "y_plot")
}

export function AnyPlotly({
  data, onChange, value
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  return VisualizationUI({data, onChange, value}, "")
}

function VisualizationUI({
  data, onChange, ...props
}: WithInputParams<{}, DataEnvelope<any>, any>, setPlotType: string = '') {
  const preset = useMemo(() => data && data['preset'], [data]);
  const hide = useMemo(() => preset && Object.keys(preset), [preset]);
  const defaultValue = whichDefaults(setPlotType, preset);
  const value = useSetsDefault(defaultValue, props.value, onChange);
  const [showAdvanced, setShowAdvanced] = useState(false);
  const plotType = (value && value['plot_type']) ? value['plot_type'] as string : '';

  const data_frame: DataEnvelope<any> = useMemo(() => {
    if (data == null) return {};
    if (data['data_frame'] == null) return {};
    return data['data_frame']
  }, [data]);
  
  const columns: string[] = useMemo(() => {
    if (data_frame == null || data_frame == {}) return [];
    return Object.keys(data_frame)
  }, [data_frame]);
  
  const updatePlotType = (newType: string, key: string) => {
    onChange(some(whichDefaults(newType, preset)))
  }

  const updateValue = (newValue: any, key: string, prevValues = {...value}) => {
    prevValues[key] = newValue;
    onChange(some(prevValues));
  };

  const extra_inputs = useExtraInputs(columns, data_frame);

  // Component set constructor
  const component_use = (key: string, value: any, extra_inputs: any) => {
    
    const comp_use: Function = comps[key]
    return(
      comp_use(key, updateValue, value, ...extra_inputs)
    )
  };

  // Plot Options
  const base = useMemo(() => {
    return (plotType!='') ? input_sets[plotType]['main'] as string[] : [] as string[]
  }, [plotType])
  
  const shownSetupValues = useMemo(() => {
    let initial = showAdvanced ? {...value} : pick(value, ...base)
    if (Object.keys(initial).includes('plot_type')) {delete initial['plot_type']}
    return remove_hidden(initial, hide);
  }, [value, showAdvanced, hide])
  
  const showHide: string = useMemo(() => {
    return (showAdvanced ? 'Hide' : 'Show')
  }, [showAdvanced])
  
  let inner = null;
  if (plotType != '') {
    
    inner = (
      <div>
        {Object.entries(shownSetupValues).map(([key, val]) => {
          return component_use(key, val, extra_inputs[key])
        })}
        <Button
          variant="contained"
          color="primary"
          onClick={() => {setShowAdvanced(!showAdvanced)}}
          >
          {showHide} Advanced Options
        </Button>
      </div>
    )
  }
  
  const pickPlot = (setPlotType!='') ? null : (
    <div>
      {dropdownInput(
      'plot_type', updatePlotType, plotType, "Plot Type:",
      Object.keys(input_sets), false)}
      <hr/>
    </div>
  )
  
  // console.log(props.value);
  
  return (
    <div>
      {pickPlot}
      {inner}
    </div>
  );

};

const remove_hidden = (vals: DataEnvelope<any>, hide: string[] | null | undefined) => {
  
  let values = {...vals};

  if (hide == null || hide.length === 0 || values.length === 0) {
    return values;
  }
  
  const keys = Object.keys(values)
  for (let ind = 0; ind < keys.length; ind++) {
    if (hide.includes(keys[ind])) delete values[keys[ind]];
  }
  return values;
};

const input_sets: DataEnvelope<DataEnvelope<string[]|DataEnvelope<any>>> = {
  'scatter_plot': {
    'main': ["x_by", "y_by", "color_by"],
    'adv': ['size', 'plot_title', 'legend_title', 'xlab', 'ylab', 'color_order', 'order_when_continuous_color', 'x_scale', 'y_scale', 'rows_use']
    //'default_adjust': {'color_by': "make"}
  },
  'bar_plot': {
    'main': ["x_by", "y_by", "scale_by"],
    'adv': ['plot_title', 'legend_title', 'xlab', 'ylab', 'rows_use']
  },
  'y_plot': {
    'main': ["x_by", "y_by", "plots"],
    'adv': ["color_by", 'plot_title', 'legend_title', 'xlab', 'ylab', 'y_scale', 'rows_use']
    //'default_adjust': {'color_by': "make"}
  }
}

const defaults: DataEnvelope<any> = {
  'x_by': null,
  'y_by': null,
  'plots': ['box', 'violin'],
  'color_by': 'make',
  'scale_by': 'fraction',
  'size': 5,
  'plot_title': 'make',
  'legend_title': 'make',
  'xlab': 'make',
  'ylab': 'make',
  'color_order': 'increasing',
  'order_when_continuous_color': false,
  'x_scale': 'as is',
  'y_scale': 'as is',
  'rows_use': {}
};

function whichDefaults(plotType: string, preset: DataEnvelope<any> | null | undefined) {
  if (plotType == '') return {plot_type: plotType}
    
  const inputs = input_sets[plotType]['main'].concat(input_sets[plotType]['adv'])

  let initial_vals = {...defaults};
  
  // Remove input:value pairs that aren't in this Viz type
  const def_keys = Object.keys(defaults);
  for (let ind = 0; ind < def_keys.length; ind++) {
    if (!inputs.includes(def_keys[ind])) delete initial_vals[def_keys[ind]];
  }
  
  // Replace any values if different default given for this plot type
  if (input_sets[plotType]['default_adjust'] != null) {
    const new_def_keys = Object.keys(input_sets[plotType]['default_adjust'])
    const new_def_vals = Object.values(input_sets[plotType]['default_adjust'])
    for (let ind = 0; ind < new_def_keys.length; ind++) {
      initial_vals[new_def_keys[ind]]=new_def_vals[ind];
    }
  }

  // Add preset values / replace any default values.
  if (preset != null) {
    const pre_keys = Object.keys(preset);
    for (let ind = 0; ind < pre_keys.length; ind++) {
      initial_vals[pre_keys[ind]] = preset[pre_keys[ind]];
    }
  }

  return {plot_type: plotType, ...initial_vals};
}

function useExtraInputs(options: string[], full_data: DataEnvelope<any>) {
  const extra_inputs: DataEnvelope<any[]> = useMemo(() => {
    return {
      // label, then for any extras
      'plot_title': ['Plot Title'],
      'legend_title': ['Legend Title'],
      'xlab': ['X-Axis Title'],
      'ylab': ['Y-Axis Title'],
      'x_by': ['X-Axis Data', options, false],
      'y_by': ['Y-Axis Data', options, false],
      'color_by': ['Color Data', ['make'].concat(options), false],
      'plots': ['Data Representations', ['violin', 'box']],
      'color_order': ['Point render order', ['increasing', 'decreasing', 'unordered']],
      'order_when_continuous_color': ['Follow selected render ordering when color is continuous?'],
      'size': ['Point Size', 0.1, 50],
      'scale_by': ['Scale Y by counts or fraction?', ['counts', 'fraction']],
      'x_scale': ['Adjust scaling of the X-Axis?', ['as is', 'log10', 'log10(val+1)']],
      'y_scale': ['Adjust scaling of the Y-Axis?', ['as is', 'log10', 'log10(val+1)']],
      'rows_use': ['Focus on a subset of the incoming data?', full_data, false, "secondary"]
    }
  }, [options, full_data]);

  return extra_inputs;
}

const comps: DataEnvelope<Function> = {
  'plot_title': stringPiece,
  'legend_title': stringPiece,
  'xlab': stringPiece,
  'ylab': stringPiece,
  'x_by': dropdownPiece,
  'y_by': dropdownPiece,
  'color_by': dropdownPiece,
  'plots': MultiselectPiece,
  'color_order': dropdownPiece,
  'order_when_continuous_color': checkboxPiece,
  'size': sliderPiece,
  'scale_by': dropdownPiece,
  'x_scale': dropdownPiece,
  'y_scale': dropdownPiece,
  'rows_use': subsetDataFramePiece
}