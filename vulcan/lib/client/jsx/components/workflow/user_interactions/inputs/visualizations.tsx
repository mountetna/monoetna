// Input component that takes nested object
//   and shows the keys one level at a time.
// Returns the last "Leaf" that the user selects.
import React, {Dispatch, PropsWithChildren, useMemo, useState} from 'react';
import * as _ from 'lodash';

import {DataEnvelope, WithInputParams} from './input_types';
import { useSetsDefault } from './useSetsDefault';
import { some } from '../../../../selectors/maybe';
import { Accordion, AccordionDetails, AccordionSummary, Grid, Typography } from '@material-ui/core';
import { pick } from 'lodash';
import { key_wrap, stringPiece, dropdownPiece, multiselectPiece, checkboxPiece, sliderPiece } from './user_input_pieces';
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
  'continuous_cols': A vector of column names of data_frame which contain continuous data.
  'discrete_cols': A vector of column names of data_frame which contain continuous data.
  (optional) 'preset': A hash where keys = input_names that should be hidden from the user & values = the preset value that should be used for that input. E.g. The x_by and y_by inputs are hardset within the umap workflow as '0' and '1', so a user has no choice here and should not be able to adjust those fields!
Note: 'data_frame', 'continuous_cols', and 'discrete_cols' can all be generated in the upstream python script by passing the dataframe through 'output_dataframe_and_types()' of from archimedes.functions.plotting.
  
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
  return VisualizationUI({data, onChange, value}, null)
}

function VisualizationUI({
  data, onChange, ...props
}: WithInputParams<{}, DataEnvelope<any>, any>, setPlotType: string | null = null) {
  const preset = useMemo(() => data && data['preset'], [data]);
  const hide = useMemo(() => preset && Object.keys(preset), [preset]);
  const defaultValue = whichDefaults(setPlotType, preset);
  const value = useSetsDefault(defaultValue, props.value, onChange);
  const [expandedDrawers, setExpandedDrawers] = useState(['primary features']);
  const plotType = (value && value['plot_type']) ? value['plot_type'] as string : null;

  const data_frame: DataEnvelope<any> = useMemo(() => {
    if (data == null) return {};
    if (data['data_frame'] == null) return {};
    return data['data_frame']
  }, [data]);
  
  const columns: string[] = useMemo(() => {
    if (data_frame == null || data_frame == {}) return [];
    return Object.keys(data_frame)
  }, [data_frame]);
  
  const continuous_columns: string[] = useMemo(() => {
    if (data == null) return [];
    if (data['continuous_cols'] == null) return columns; // Should build a warning here instead?
    return data['continuous_cols']
  }, [data]);
  const discrete_columns: string[] = useMemo(() => {
    if (data == null) return [];
    if (data['discrete_cols'] == null) return columns; // Should build a warning here instead?
    return data['discrete_cols']
  }, [data]);

  const extra_inputs = useExtraInputs(columns, data_frame, plotType, continuous_columns, discrete_columns);
  
  const shownSetupValues = useMemo(() => {
    if (plotType==null) return {}
    const initial = {...value}
    if (Object.keys(initial).includes('plot_type')) {delete initial['plot_type']}
    return remove_hidden(initial, hide);
  }, [value, hide])
  
  const updatePlotType = (newType: string, key: string) => {
    onChange(some(whichDefaults(newType, preset)))
  }

  const updateValue = (newValue: any, key: string, prevValues = {...value}) => {
    prevValues[key] = newValue;
    onChange(some(prevValues));
  };
  
  function toggleDrawerExpansion(drawerTitle: string) {
    if (expandedDrawers.includes(drawerTitle)) {
      setExpandedDrawers(expandedDrawers.filter(t => t!=drawerTitle))
    } else {
      setExpandedDrawers(expandedDrawers.concat(drawerTitle))
    }
  }

  // Components
  const pickPlot = (setPlotType!=null) ? null : (
    <div>
      {dropdownPiece(
      'plot_type', updatePlotType, plotType, "Plot Type",
      Object.keys(input_sets), false)}
    </div>
  )
  
  const inner = (plotType == null) ? null : (
    <Grid container direction='column'>
      {Object.entries(input_sets[plotType]).map(([group_name, val]) => {
        const group_values: DataEnvelope<any> = pick(shownSetupValues,val)
        const open = expandedDrawers.includes(group_name);
        // console.log('group_values', group_values)
        return (Object.keys(group_values).length > 0) ? <InputWrapper key={group_name} title={group_name} values={group_values} open={open} toggleOpen={toggleDrawerExpansion}>
          {Object.entries(group_values).map(([key, val]) => {
            return <ComponentUse key={key} k={key} value={val} extra_inputs={extra_inputs[key]} updateValue={updateValue}/>
          })}
          </InputWrapper> : null
      })}
    </Grid>
  )
  
  // console.log(props.value);
  
  return (
    <div key='VizUI'>
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

const input_sets: DataEnvelope<DataEnvelope<string[]>> = {
  'scatter_plot': {
    'primary features': ["x_by", "y_by", "color_by", 'size'],
    'titles': ['plot_title', 'legend_title', 'xlab', 'ylab'],
    'point rendering': ['color_order', 'order_when_continuous_color'],
    'coordinates': ['x_scale', 'y_scale'],
    'data focus': ['rows_use']
    //'default_adjust': {'color_by': "make"}
  },
  'bar_plot': {
    'primary features': ["x_by", "y_by", "scale_by"],
    'titles': ['plot_title', 'legend_title', 'xlab', 'ylab'],
    'data focus': ['rows_use']
  },
  'y_plot': {
    'primary features': ["x_by", "y_by", "plots", "color_by"],
    'titles': ['plot_title', 'legend_title', 'xlab', 'ylab'],
    'coordinates': ['y_scale'],
    'data focus': ['rows_use']
    //'default_adjust': {'color_by': "make"}
  }
}

const input_constraints: DataEnvelope<DataEnvelope<"continuous"|"discrete">> = {
  'scatter_plot': {
    'x_by': "continuous",
    'y_by': "continuous"
  },
  'bar_plot': {
    'x_by': "discrete",
    'y_by': "discrete"
  },
  'y_plot': {
    'x_by': "discrete",
    'y_by': "continuous"
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
  'x_scale': 'linear',
  'y_scale': 'linear',
  'rows_use': {}
};

function whichDefaults(plotType: string|null, preset: DataEnvelope<any> | null | undefined) {
  if (plotType == null) return {plot_type: plotType}
  
  const inputs = Object.values(input_sets[plotType]).flat()

  let initial_vals = {...defaults};
  
  // Remove input:value pairs that aren't in this Viz type
  const def_keys = Object.keys(defaults);
  for (let ind = 0; ind < def_keys.length; ind++) {
    if (!inputs.includes(def_keys[ind])) delete initial_vals[def_keys[ind]];
  }
  
  // // Replace any values if different default given for this plot type
  // Broken during a redesign, but should be restorable easily!
  // if (input_sets[plotType]['default_adjust'] != null) {
  //   const new_def_keys = Object.keys(input_sets[plotType]['default_adjust'])
  //   const new_def_vals = Object.values(input_sets[plotType]['default_adjust'])
  //   for (let ind = 0; ind < new_def_keys.length; ind++) {
  //     initial_vals[new_def_keys[ind]]=new_def_vals[ind];
  //   }
  // }

  // Add preset values / replace any default values.
  if (preset != null) {
    const pre_keys = Object.keys(preset);
    for (let ind = 0; ind < pre_keys.length; ind++) {
      initial_vals[pre_keys[ind]] = preset[pre_keys[ind]];
    }
  }

  return {plot_type: plotType, ...initial_vals};
}

function useExtraInputs(
  options: string[], full_data: DataEnvelope<any>,
  plot_type: string | null,
  continuous: string[], discrete: string[],
  constraints: DataEnvelope<DataEnvelope<"continuous"|"discrete">> = input_constraints
  ) {
  
  function get_options(input_name: string) {
    if (plot_type==null) return options
    if (Object.keys(constraints[plot_type]).includes(input_name)) {
      if (constraints[plot_type][input_name]=="continuous") return continuous
      if (constraints[plot_type][input_name]=="discrete") return discrete
    }
    return options
  }
  
  const extra_inputs: DataEnvelope<any[]> = useMemo(() => {
    return {
      // label, then for any extras
      'plot_title': ['Plot Title'],
      'legend_title': ['Legend Title'],
      'xlab': ['X-Axis Title'],
      'ylab': ['Y-Axis Title'],
      'x_by': ['X-Axis Data', get_options('x_by'), false],
      'y_by': ['Y-Axis Data', get_options('y_by'), false],
      'color_by': ['Color Data', ['make'].concat(get_options('color_by')), false],
      'plots': ['Data Representations', ['violin', 'box']],
      'color_order': ['Point render order', ['increasing', 'decreasing', 'unordered']],
      'order_when_continuous_color': ['Follow selected render ordering when color is continuous?'],
      'size': ['Point Size', 0.1, 50],
      'scale_by': ['Scale Y by counts or fraction', ['counts', 'fraction']],
      'x_scale': ['Adjust scaling of the X-Axis', ['linear', 'log10', 'log10(val+1)']],
      'y_scale': ['Adjust scaling of the Y-Axis', ['linear', 'log10', 'log10(val+1)']],
      'rows_use': ['Focus on a subset of the incoming data', full_data, false, "secondary"]
    }
  }, [options, plot_type, constraints, continuous, discrete]);

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
  'plots': multiselectPiece,
  'color_order': dropdownPiece,
  'order_when_continuous_color': checkboxPiece,
  'size': sliderPiece,
  'scale_by': dropdownPiece,
  'x_scale': dropdownPiece,
  'y_scale': dropdownPiece,
  'rows_use': subsetDataFramePiece
}

function InputWrapper({title, values, open, toggleOpen, children}: PropsWithChildren<{title:string, values: DataEnvelope<any>, open: boolean, toggleOpen: Dispatch<string>}>) {
  return (
    <Accordion elevation={0} expanded={open} onChange={ () => toggleOpen(title) }>
      <AccordionSummary
        style={{background: '#eee', cursor: 'pointer', minHeight: '32px',  height: '32px'}}>
        <Typography>{title}</Typography>
      </AccordionSummary>
      <AccordionDetails
        style={{padding: 0, paddingLeft: 15, borderBottom: '1px solid #eee'}}>
        <Grid container direction='column'>
          {children}
        </Grid>
      </AccordionDetails>
    </Accordion>
  )
}

const ComponentUse = ({k, value, extra_inputs, updateValue}: {k: string, value: any, extra_inputs: any, updateValue: Function}) => {
    
  const comp_use: Function = comps[k]
  return(
    comp_use(k, updateValue, value, ...extra_inputs)
  )
};