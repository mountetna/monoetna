// Input component that takes nested object
//   and shows the keys one level at a time.
// Returns the last "Leaf" that the user selects.
import React, {Dispatch, PropsWithChildren, useMemo, useState} from 'react';
import * as _ from 'lodash';

import {DataEnvelope, WithInputParams} from './input_types';
import {useSetsDefault} from './useSetsDefault';
import {some} from '../../../../selectors/maybe';
import {
  Accordion,
  AccordionDetails,
  AccordionSummary,
  Grid,
  Typography
} from '@material-ui/core';
import {pick} from 'lodash';
import {
  key_wrap,
  stringPiece,
  dropdownPiece,
  multiselectPiece,
  checkboxPiece,
  sliderPiece,
  reductionSetupPiece,
  nestedDropdownPiece
} from './user_input_pieces';
import {subsetDataFramePiece} from './subsetDataFrame_piece';
import {ReorderPiece} from './reorder_piece';

/*
Docmentation last updated: Apr 15, 2022

This UI is closely tied to archimedes/functions/plotting/*_plotly functions.

CWL Call: (parenthesis = optional)

  <step(s)_before must have>
      out: [data_frame, continuous_cols, discrete_cols, (preset)]

  user_plot_setup:
    run: ui-queries/any-viz.cwl (Alternately, to lock in a plot type: 'ui-queries/scatter-plotly.cwl', 'ui-queries/y-plotly.cwl', or 'ui-queries/bar-plotly.cwl')
    label: 'Set plot options'
    doc: "Selections here adjust how the plot will be generated. For addtional details, see the 'Visualization with Vulcan' section of the Vulcan's documentation, acccessible via the 'Help' button at the top of this page. (<if using 'presets' or a set plot type, please describe! Example, scRNAseq.cwl umap:> This particular instance of the Plot Configuration Interface constitutes a version with preset values for plot-type (scatter_plot), X-Axis Data (UMAP_1), Y-Axis Data (UMAP_2), and Color Data (chosen above).)"
    in:
      data_frame: <in_step>/data_frame
      continuous_cols: <in_step>/continuous_cols
      discrete_cols: <in_step>/discrete_cols
      (optional)
      preset: <in_step>/preset
    out: [plot_setup]
  make_umap:
    run: scripts/make_plot.cwl
    label: 'Create UMAP plot'
    in:
      plot_setup: user_plot_setup/plot_setup
      data_frame: <in_step>/data_frame
    out: [plot.json, plot.png]

Creation from python:
  Imports:
    from archimedes.functions.dataflow import output_path
    from archimedes.functions.plotting import output_column_types
  dataframe:
    'df.to_json(output_path("data_frame"))', where df is a pandas.DataFrame with attributes in columns and data points/obsevations in rows.
  continuous_cols, discrete_cols:
    'output_column_types(df, "continuous_cols", "discrete_cols")', with the same df used for 'dataframe'
  preset:
    'output_json(preset, "preset")' where preset is a dict which one could pass to the target plotter function (in 'archimedes/archimedes/functions/plotting') via 'viz_fxn(**preset)'
    Example, a umap where color_by is pre-selected in an upstream ui, and umap embeddings are stored in columns named '0' and '1':
      color_by = input_var('color_by')
      plot_preset = {
        'x_by': '0',
        'y_by': '1',
        'color_by': color_by,
        'xlab': 'UMAP_1',
        'ylab': 'UMAP_2'
      }
      output_json(preset, "preset")

JSX:
  'data', object with slots:
    data_frame:
      A hash representing the data_frame that will be used to make a plot. keys = column names, values = data points.
    continuous_cols:
      A vector of column names of data_frame which contain continuous data.
    discrete_cols:
      A vector of column names of data_frame which contain discrete data.
    ?preset:
      A hash where keys = input pieces that should be hidden from the user & values = the preset value that should be used for that input. E.g. The x_by and y_by inputs are hardset within the umap workflow as '0' and '1', so a user has no choice here and should not be able to adjust those fields!

  'value', object where (key,val) pairs mostly equate to (key,val) pairs that could be splatted into an archimedes visualization funciton.

  Major design notes:
    - Organized around having pairs of input-names and associated component-setups.
    - An 'input_sets' object defines the set of inputs that are used for a given visualization type.  Thus, after a plot-type choice, necessary inputs are known and the set of necessary ui pieces can be compiled.
    - (key,val) pairs of the optional 'presets' input will cause any inputs' component-setup to be removed from the displayed set, while also providing the value to give to 'value[key]'.

*/

type optionSet = string[]
type nestedOptionSet = DataEnvelope<DataEnvelope<DataEnvelope<null>|null>|null>

const defaults_plotly: DataEnvelope<any> = {
  x_by: null,
  y_by: null,
  plots: ['box', 'violin'],
  color_by: 'make',
  scale_by: 'fraction',
  size: 5,
  plot_title: 'make',
  legend_title: 'make',
  xlab: 'make',
  ylab: 'make',
  color_order: 'increasing',
  order_when_continuous_color: false,
  x_scale: 'linear',
  y_scale: 'linear',
  rows_use: {},
  x_order: 'increasing',
  y_order: 'increasing'
};

const defaults_dittoseq: DataEnvelope<any> = {
  x_by: null,
  y_by: null,
  var: null,
  group_by: null,
  color_by: null,
  plots: ['violin', 'jitter'],
  scale_by: 'fraction',
  size: 1,
  plot_title: 'make',
  legend_title: 'make',
  xlab: 'make',
  ylab: 'make',
  color_order: 'increasing',
  x_scale: 'linear',
  y_scale: 'linear',
  cells_use: {},
  reduction_setup: ['_Recommended_', '1', '2'],
  do_hover: true
};

const remove_hidden = (
  vals: DataEnvelope<any>,
  hide: string[] | null | undefined
) => {
  let values = {...vals};

  if (hide == null || hide.length === 0 || values.length === 0) {
    return values;
  }

  const keys = Object.keys(values);
  for (let ind = 0; ind < keys.length; ind++) {
    if (hide.includes(keys[ind])) delete values[keys[ind]];
  }
  return values;
};

const input_sets_plotly: DataEnvelope<DataEnvelope<string[]>> = {
  scatter_plot: {
    'primary features': ['x_by', 'y_by', 'color_by', 'size'],
    titles: ['plot_title', 'legend_title', 'xlab', 'ylab'],
    coordinates: ['x_scale', 'y_scale'],
    'data focus': ['rows_use', 'color_order', 'order_when_continuous_color']
    //'default_adjust': {'color_by': "make"}
  },
  bar_plot: {
    'primary features': ['x_by', 'y_by', 'scale_by'],
    titles: ['plot_title', 'legend_title', 'xlab', 'ylab'],
    'data focus': ['rows_use']
  },
  y_plot: {
    'primary features': ['x_by', 'y_by', 'plots', 'color_by'],
    titles: ['plot_title', 'legend_title', 'xlab', 'ylab'],
    coordinates: ['y_scale'],
    'data focus': ['rows_use', 'x_order']
    //'default_adjust': {'color_by': "make"}
  }
};

const input_sets_dittoseq: DataEnvelope<DataEnvelope<string[]>> = {
  dittoDimPlot: {
    'primary features': ['color_by', 'size', 'reduction_setup'],
    titles: ['plot_title', 'legend_title', 'xlab', 'ylab'],
    'data focus': ['color_order', 'cells_use'],
    'output style': ['do_hover']
  },
  dittoScatterPlot: {
    'primary features': ['x_by', 'y_by', 'color_by', 'size'],
    titles: ['plot_title', 'legend_title', 'xlab', 'ylab'],
    'data focus': ['color_order', 'cells_use'],
    'output style': ['do_hover']
  },
  dittoBarPlot: {
    'primary features': ['var', 'group_by', 'scale_by'],
    titles: ['plot_title', 'legend_title', 'xlab', 'ylab'],
    'output style': ['do_hover'],
    'data focus': ['cells_use']
  },
  dittoPlot: {
    'primary features': ['var', 'group_by', 'plots', 'color_by'],
    titles: ['plot_title', 'legend_title', 'xlab', 'ylab'],
    'data focus': ['cells_use']
  }
};

function whichDefaults(
  plotType: string | null,
  preset: DataEnvelope<any> | null | undefined,
  defaults: DataEnvelope<any>,
  input_sets: DataEnvelope<DataEnvelope<string[]>>
) {
  if (plotType == null) return {plot_type: plotType};

  const inputs = Object.values(input_sets[plotType]).flat();

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

const input_constraints: DataEnvelope<DataEnvelope<'continuous' | 'discrete'>> =
  {
    scatter_plot: {
      x_by: 'continuous',
      y_by: 'continuous'
    },
    bar_plot: {
      x_by: 'discrete',
      y_by: 'discrete'
    },
    y_plot: {
      x_by: 'discrete',
      y_by: 'continuous',
      color_by: 'discrete'
    },
    dittoDimPlot: {
      split_by: 'discrete'
    },
    dittoScatterPlot: {
      x_by: 'continuous',
      y_by: 'continuous',
      split_by: 'discrete'
    },
    dittoBarPlot: {
      group_by: 'discrete',
      var: 'discrete',
      split_by: 'discrete'
    },
    dittoPlot: {
      group_by: 'discrete',
      var: 'continuous',
      color_by: 'discrete',
      split_by: 'discrete'
    }
  };

function useExtraInputs(
  options: string[] | nestedOptionSet,
  full_data: DataEnvelope<any>,
  plot_type: string | null,
  continuous: string[] | nestedOptionSet,
  discrete: string[] | nestedOptionSet,
  x_by: string | null,
  y_by: string | null,
  color_by: string | null,
  reduction_opts: DataEnvelope<number[]> | null = null,
  constraints: DataEnvelope<
    DataEnvelope<'continuous' | 'discrete'>
  > = input_constraints
) {
  function get_options(input_name: string) {
    if (plot_type == null) return options;
    if (Object.keys(constraints[plot_type]).includes(input_name)) {
      if (constraints[plot_type][input_name] == 'continuous') return continuous;
      if (constraints[plot_type][input_name] == 'discrete') return discrete;
    }
    return options;
  }

  const extra_inputs: DataEnvelope<any[]> = useMemo(() => {
    function is_ditto() {
      return plot_type!=null && plot_type.includes('ditto')
    }
    function add_make(add_to: string[] | nestedOptionSet, for_plot_types: string[]) {
      if (plot_type==null || !for_plot_types.includes(plot_type)) return add_to
      if (Array.isArray(add_to)) return ['make'].concat(add_to)
      let output = {...add_to}
      output['make']=null
      return output
    }
    return {
      // label, then for any extras
      plot_title: ['Plot Title'],
      legend_title: ['Legend Title'],
      xlab: ['X-Axis Title'],
      ylab: ['Y-Axis Title'],
      x_by: ['X-Axis Data', get_options('x_by'), false],
      y_by: ['Y-Axis Data', get_options('y_by'), false],
      var: ['Primary Data', get_options('var'), false],
      group_by: ['Groupings Data (often the x-axis)', get_options('group_by'), false],
      color_by: ['Color Data', add_make(get_options('color_by'), ['scatter_plot', 'y_plot', 'dittoPlot']), false],
      plots: ['Data Representations', ['violin', 'box']],
      color_order: [
        'Point Render & (discrete) Color Assignment Order',
        full_data,
        color_by,
        discrete,
        is_ditto()
      ],
      order_when_continuous_color: [
        'Follow selected render ordering when color is continuous?'
      ],
      size: !is_ditto ? ['Point Size', 0.1, 50, undefined] : ['Point Size', 0.1, 25, 0.1],
      scale_by: [
        'Scale Y by counts or fraction',
        ['counts', 'fraction'],
        true,
        200,
        x_by == y_by
      ],
      x_scale: [
        'Adjust scaling of the X-Axis',
        ['linear', 'log10', 'log10(val+1)']
      ],
      y_scale: [
        'Adjust scaling of the Y-Axis',
        ['linear', 'log10', 'log10(val+1)']
      ],
      rows_use: [
        'Focus on a subset of the incoming data',
        full_data,
        false,
        'secondary'
      ],
      cells_use: [
        'Focus on a subset of cells',
        full_data,
        false,
        'secondary',
        continuous,
        'Discrete Cell_Metadata'
      ],
      x_order: ['Order of X-Axis Groupings', full_data, x_by, discrete],
      y_order: ['Order of Y-Axis Groupings', full_data, y_by, discrete],
      reduction_setup: [
        ['Dimensionality Reduction (DR)', 'x-axis DR Compenent #', 'y-axis DR Component #'],
        reduction_opts
      ],
      do_hover: ['Output as interactive (not flat) image?']
    };
  }, [
    options,
    plot_type,
    constraints,
    continuous,
    discrete,
    x_by,
    y_by,
    color_by,
    reduction_opts
  ]);

  return extra_inputs;
}

function InputWrapper({
  title,
  values,
  open,
  toggleOpen,
  children
}: PropsWithChildren<{
  title: string;
  values: DataEnvelope<any>;
  open: boolean;
  toggleOpen: Dispatch<string>;
}>) {
  return (
    <Accordion elevation={0} expanded={open} onChange={() => toggleOpen(title)}>
      <AccordionSummary
        style={{
          background: '#eee',
          cursor: 'pointer',
          minHeight: '32px',
          height: '32px'
        }}
      >
        <Typography>{title}</Typography>
      </AccordionSummary>
      <AccordionDetails
        style={{padding: 0, paddingLeft: 15, borderBottom: '1px solid #eee'}}
      >
        <Grid container direction='column'>
          {children}
        </Grid>
      </AccordionDetails>
    </Accordion>
  );
}

const components_plotly: DataEnvelope<Function> = {
  plot_title: stringPiece,
  legend_title: stringPiece,
  xlab: stringPiece,
  ylab: stringPiece,
  x_by: dropdownPiece,
  y_by: dropdownPiece,
  color_by: dropdownPiece,
  plots: multiselectPiece,
  color_order: ReorderPiece,
  order_when_continuous_color: checkboxPiece,
  size: sliderPiece,
  scale_by: dropdownPiece,
  x_scale: dropdownPiece,
  y_scale: dropdownPiece,
  rows_use: subsetDataFramePiece,
  x_order: ReorderPiece,
  y_order: ReorderPiece
};

const components_dittoseq: DataEnvelope<Function> = {
  plot_title: stringPiece,
  legend_title: stringPiece,
  xlab: stringPiece,
  ylab: stringPiece,
  x_by: nestedDropdownPiece,
  y_by: nestedDropdownPiece,
  var: nestedDropdownPiece,
  group_by: nestedDropdownPiece,
  color_by: nestedDropdownPiece,
  plots: multiselectPiece,
  color_order: ReorderPiece,
  size: sliderPiece,
  scale_by: dropdownPiece,
  x_scale: dropdownPiece,
  y_scale: dropdownPiece,
  cells_use: subsetDataFramePiece,
  x_order: ReorderPiece,
  y_order: ReorderPiece,
  reduction_setup: reductionSetupPiece,
  do_hover: checkboxPiece,
};

const ComponentUse = ({
  k,
  value,
  extra_inputs,
  updateValue,
  comps
}: {
  k: string;
  value: any;
  extra_inputs: any;
  updateValue: Function;
  comps: DataEnvelope<Function>;
}) => {
  const comp_use: Function = comps[k];
  return comp_use(k, updateValue, value, ...extra_inputs);
};

function VisualizationUI(
  {data, onChange, ...props}: WithInputParams<{}, DataEnvelope<any>, any>,
  setPlotType: string | null = null,
  defaults: DataEnvelope<any>,
  input_sets: DataEnvelope<DataEnvelope<string[]>>,
  components: DataEnvelope<Function>
) {
  const preset = useMemo(() => data && data['preset'], [data]);
  const hide = useMemo(() => preset && Object.keys(preset), [preset]);
  const defaultValue = whichDefaults(setPlotType, preset, defaults, input_sets);
  const value = useSetsDefault(defaultValue, props.value, onChange);
  const [expandedDrawers, setExpandedDrawers] = useState(['primary features']);
  const plotType =
    value && value['plot_type'] ? (value['plot_type'] as string) : null;

  const data_frame: DataEnvelope<any> = useMemo(() => {
    if (data == null) return {};
    if (data['data_frame'] == null) return {};
    return data['data_frame'];
  }, [data]);

  const df_columns: string[] = useMemo(() => {
    if (data_frame == null || 0 === Object.keys(data_frame).length) return [];
    return Object.keys(data_frame);
  }, [data_frame]);

  const continuous_columns: string[] | nestedOptionSet = useMemo(() => { // NOTE: these data don't necessarily need to be contained within the given data_frame (to accomodate for visualizing genomics data)
    if (data == null) return [];
    if (data['continuous_cols'] == null) return df_columns; // Should build a warning here instead?
    return data['continuous_cols'];
  }, [data]);
  const discrete_columns: string[] | nestedOptionSet = useMemo(() => {
    if (data == null) return [];
    if (data['discrete_cols'] == null) return df_columns; // Should build a warning here instead?
    return data['discrete_cols'];
  }, [data]);

  const columns: string[] | nestedOptionSet = (data != null && data['all_cols'] != null) ? data['all_cols'] : df_columns

  const reduction_opts: DataEnvelope<string[]> | null = useMemo(() => {
    if (data == null) return null;
    if (data['reduction_opts'] == null) return null;
    return data['reduction_opts'];
  }, [data])
  // console.log({reduction_opts});

  const x_by =
    value && Object.keys(value).includes('x_by')
      ? (value.x_by as string | null)
      : null;
  const y_by =
    value && Object.keys(value).includes('y_by')
      ? (value.y_by as string | null)
      : null;
  const color_by =
    value && Object.keys(value).includes('color_by')
      ? (value.color_by as string | null)
      : null;
  const extra_inputs = useExtraInputs(
    columns,
    data_frame,
    plotType,
    continuous_columns,
    discrete_columns,
    x_by,
    y_by,
    color_by,
    reduction_opts
  );

  const shownSetupValues = useMemo(() => {
    if (plotType == null) return {};
    const initial = {...value};
    if (Object.keys(initial).includes('plot_type')) {
      delete initial['plot_type'];
    }
    return remove_hidden(initial, hide);
  }, [value, hide]);

  const updatePlotType = (newType: string, key: string) => {
    onChange(some(whichDefaults(newType, preset, defaults, input_sets)));
  };

  const updateValue = (newValue: any, key: string, prevValues = {...value}) => {
    prevValues[key] = newValue;
    onChange(some(prevValues));
  };

  function toggleDrawerExpansion(drawerTitle: string) {
    if (expandedDrawers.includes(drawerTitle)) {
      setExpandedDrawers(expandedDrawers.filter((t) => t != drawerTitle));
    } else {
      setExpandedDrawers(expandedDrawers.concat(drawerTitle));
    }
  }

  // Components
  const pickPlot =
    setPlotType != null ? null : (
      <div>
        {dropdownPiece(
          'plot_type',
          updatePlotType,
          plotType,
          'Plot Type',
          Object.keys(input_sets),
          false
        )}
      </div>
    );

  const inner =
    plotType == null ? null : (
      <Grid container direction='column'>
        {Object.entries(input_sets[plotType]).map(([group_name, val]) => {
          const group_values: DataEnvelope<any> = pick(shownSetupValues, val);
          const open = expandedDrawers.includes(group_name);
          // console.log('group_values', group_values)
          return Object.keys(group_values).length > 0 ? (
            <InputWrapper
              key={group_name}
              title={group_name}
              values={group_values}
              open={open}
              toggleOpen={toggleDrawerExpansion}
            >
              {Object.entries(group_values).map(([key, val]) => {
                return (
                  <ComponentUse
                    key={key}
                    k={key}
                    value={val}
                    extra_inputs={extra_inputs[key]}
                    updateValue={updateValue}
                    comps={components}
                  />
                );
              })}
            </InputWrapper>
          ) : null;
        })}
      </Grid>
    );

  // console.log(props.value);

  return (
    <div key='VizUI'>
      {pickPlot}
      {inner}
    </div>
  );
}

export function ScatterPlotly({
  data,
  onChange,
  value
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  return VisualizationUI({data, onChange, value}, 'scatter_plot', defaults_plotly, input_sets_plotly, components_plotly);
}

export function BarPlotly({
  data,
  onChange,
  value
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  return VisualizationUI({data, onChange, value}, 'bar_plot', defaults_plotly, input_sets_plotly, components_plotly);
}

export function YPlotly({
  data,
  onChange,
  value
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  return VisualizationUI({data, onChange, value}, 'y_plot', defaults_plotly, input_sets_plotly, components_plotly);
}

export function AnyPlotly({
  data,
  onChange,
  value
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  return VisualizationUI({data, onChange, value}, null, defaults_plotly, input_sets_plotly, components_plotly);
}

export function DittoDimPlot({
  data,
  onChange,
  value
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  return VisualizationUI({data, onChange, value}, 'dittoDimPlot', defaults_dittoseq, input_sets_dittoseq, components_dittoseq);
}

export function DittoScatterPlot({
  data,
  onChange,
  value
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  return VisualizationUI({data, onChange, value}, 'dittoScatterPlot', defaults_dittoseq, input_sets_dittoseq, components_dittoseq);
}

export function DittoBarPlot({
  data,
  onChange,
  value
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  return VisualizationUI({data, onChange, value}, 'dittoBarPlot', defaults_dittoseq, input_sets_dittoseq, components_dittoseq);
}

export function DittoPlot({
  data,
  onChange,
  value
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  return VisualizationUI({data, onChange, value}, 'dittoPlot', defaults_dittoseq, input_sets_dittoseq, components_dittoseq);
}

export function AnyDittoSeq({
  data,
  onChange,
  value
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  return VisualizationUI({data, onChange, value}, null, defaults_dittoseq, input_sets_dittoseq, components_dittoseq);
}
