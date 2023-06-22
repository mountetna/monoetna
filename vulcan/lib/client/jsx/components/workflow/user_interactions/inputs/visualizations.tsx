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
  stringPiece,
  dropdownPiece,
  MultiselectPiece,
  checkboxPiece,
  sliderPiece,
  reductionSetupPiece,
  nestedDropdownPiece,
  MultiselectAfterDataChoicePiece_forDitto
} from './user_input_pieces';
import {subsetDataFramePiece} from './subsetDataFrame_piece';
import {ReorderCustomOnlyPiece, ReorderPiece} from './reorder_piece';

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
export type nestedOptionSet = DataEnvelope<DataEnvelope<DataEnvelope<null>|null>|null>

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
    //'default_adjust': {'color_by': 'make'}
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
    //'default_adjust': {'color_by': 'make'}
  }
};

const plot_relabels_plotly: DataEnvelope<string> = {
  'scatter_plot: Scatter Plot': 'scatter_plot',
  'bar_plot: Compositional Stacked Bar Plot': 'bar_plot',
  'y_plot: Violin and/or Box Plot': 'y_plot'
}

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
  x_order: 'increasing'
};

const redefaults_plotly: DataEnvelope<DataEnvelope<any>> = {
  scatter_plot: {
    'color_by': 'make'
  },
  y_plot: {
    'color_by': 'make'
  }
}

const input_sets_dittoseq: DataEnvelope<DataEnvelope<string[]>> = {
  dittoDimPlot: {
    'primary features': ['color_by', 'size', 'reduction_setup'],
    titles: ['plot_title', 'plot_subtitle', 'legend_title', 'xlab', 'ylab'],
    'data focus': ['color_order', 'cells_use'],
    'addition: labels': ['do_label', 'labels_highlight', 'labels_repel'],
    'other additions': ['do_contour', 'do_ellipse'],
    'output style': ['do_hover', 'legend_show']
  },
  dittoScatterPlot: {
    'primary features': ['x_by', 'y_by', 'color_by', 'size'],
    titles: ['plot_title', 'plot_subtitle', 'legend_title', 'xlab', 'ylab'],
    'data focus': ['color_order', 'cells_use'],
    'addition: labels': ['do_label', 'labels_highlight', 'labels_repel'],
    'other additions': ['do_contour', 'do_ellipse'],
    'output style': ['do_hover', 'legend_show']
  },
  dittoBarPlot: {
    'primary features': ['var', 'group_by', 'scale_by'],
    titles: ['plot_title', 'plot_subtitle', 'legend_title', 'xlab', 'ylab'],
    'output style': ['do_hover', 'legend_show'],
    'data focus': ['group_order', 'var_order', 'cells_use']
  },
  dittoPlot: {
    'primary features': ['var', 'group_by', 'plots', 'color_by'],
    titles: ['plot_title', 'plot_subtitle', 'legend_title', 'xlab', 'ylab'],
    'vlnplot tweaks': ['vlnplot_lineweight', 'vlnplot_width', 'vlnplot_scaling'],
    'boxplot tweaks': ['boxplot_width', 'boxplot_color', 'boxplot_fill', 'boxplot_lineweight'],
    'jitter tweaks': ['jitter_size', 'jitter_width', 'jitter_color'],
    'ridgeplot tweaks': ['ridgeplot_lineweight'],
    'output style': ['legend_show'],
    'data focus': ['group_order', 'cells_use'],
  },
  dittoFreqPlot: {
    'primary features': ['var', 'vars_use', 'group_by', 'sample_by', 'scale_by', 'plots', 'color_by'],
    titles: ['plot_title', 'plot_subtitle', 'legend_title', 'xlab', 'ylab'],
    'vlnplot tweaks': ['vlnplot_lineweight', 'vlnplot_width', 'vlnplot_scaling'],
    'boxplot tweaks': ['boxplot_width', 'boxplot_color', 'boxplot_fill', 'boxplot_lineweight'],
    'jitter tweaks': ['jitter_size', 'jitter_width', 'jitter_color'],
    'ridgeplot tweaks': ['ridgeplot_lineweight'],
    'output style': ['legend_show'],
    'data focus': ['group_order', 'cells_use'],
  }
};

const plot_relabels_dittoseq: DataEnvelope<string> = {
  'dittoDimPlot: UMAP, PCA, etc.': 'dittoDimPlot',
  'dittoScatterPlot: e.g. gene x gene': 'dittoScatterPlot',
  'dittoBarPlot: compositional stacked bar plot': 'dittoBarPlot',
  'dittoPlot: violin, box, or ridge plot': 'dittoPlot',
  'dittoFreqPlot: compositional violin, box, or ridge plot': 'dittoFreqPlot'
}

const defaults_dittoseq: DataEnvelope<any> = {
  x_by: null,
  y_by: null,
  var: null,
  group_by: null,
  color_by: null,
  sample_by: 'make',
  plots: ['jitter', 'vlnplot'],
  scale_by: 'fraction',
  size: 1,
  plot_title: 'make',
  plot_subtitle: 'make',
  legend_title: 'make',
  xlab: 'make',
  ylab: 'make',
  color_order: 'unordered',
  group_order: 'make',
  var_order: 'make',
  vars_use: 'make',
  x_scale: 'linear',
  y_scale: 'linear',
  cells_use: {},
  reduction_setup: ['_Recommended_', '1', '2'],
  do_hover: true,
  vlnplot_lineweight: 1,
  vlnplot_width: 1,
  vlnplot_scaling: 'area',
  boxplot_width: 0.2,
  boxplot_color: 'black',
  boxplot_fill: true,
  boxplot_lineweight: 1,
  jitter_size: 1,
  jitter_width: 0.2,
  jitter_color: 'black',
  ridgeplot_lineweight: 1,
  do_label: false,
  labels_highlight: true,
  labels_repel: true,
  do_contour: false,
  do_ellipse: false,
  legend_show: true
};

const redefaults_dittoseq: DataEnvelope<DataEnvelope<any>> = {
  dittoScatterPlot: {
    'color_by': 'make'
  },
  dittoPlot: {
    'color_by': 'make'
  },
  dittoFreqPlot: {
    'color_by': 'make',
    'plots': ['jitter', 'boxplot']
  }
}

function whichDefaults(
  plotType: string | null,
  preset: DataEnvelope<any> | null | undefined,
  defaults: DataEnvelope<any>,
  redefaults: DataEnvelope<DataEnvelope<any>>,
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

  // Replace default values with 'redefault' values.
  if (Object.keys(redefaults).includes(plotType)) {
    const new_def_keys = Object.keys(redefaults[plotType])
    const new_def_vals = Object.values(redefaults[plotType])
    for (let ind = 0; ind < new_def_keys.length; ind++) {
      initial_vals[new_def_keys[ind]]=new_def_vals[ind];
    }
  }

  // Use preset values.
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
    },
    dittoFreqPlot: {
      group_by: 'discrete',
      var: 'discrete',
      color_by: 'discrete',
      split_by: 'discrete',
      sample_by: 'discrete'
    }
  };

function useExtraInputs(
  options: string[] | nestedOptionSet,
  full_data: DataEnvelope<any>,
  plot_type: string | null,
  continuous: string[] | nestedOptionSet,
  discrete: string[],
  x_by: string | null,
  y_by: string | null,
  color_by: string | null,
  group_by: string | null,
  var_: string | null,
  reduction_opts: DataEnvelope<number[]> | null = null,
  constraints: DataEnvelope<DataEnvelope<'continuous' | 'discrete'>> = input_constraints) {
  
  function get_options(input_name: string) {
    if (plot_type == null) return options;
    if (Object.keys(constraints[plot_type]).includes(input_name)) {
      if (constraints[plot_type][input_name] == 'continuous') return continuous;
      if (constraints[plot_type][input_name] == 'discrete') return discrete;
    }
    return options;
  }
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
  
  const extra_inputs: DataEnvelope<any[]> = useMemo(() => {
    return {
      // label, then for any extras
      plot_title: ['Plot Title'],
      plot_subtitle: ['Plot Sub-Title'],
      legend_title: ['Legend Title'],
      xlab: ['X-Axis Title'],
      ylab: ['Y-Axis Title'],
      x_by: ['X-Axis Data', get_options('x_by'), false],
      y_by: ['Y-Axis Data', get_options('y_by'), false],
      var: ['Primary Data', get_options('var'), false],
      group_by: ['Groupings Data (often the x-axis)', get_options('group_by'), false],
      color_by: [
        is_ditto() ? 'Color Data (\'make\' = same as Groupings Data)' : 'Color Data',
        add_make(get_options('color_by'), ['scatter_plot', 'y_plot', 'dittoPlot', 'dittoFreqPlot', 'dittoScatterPlot']),
        false],
      sample_by: ['Sample Data (try \'_sc_seq_ids_\' if there; \'make\' = ignore)', add_make(get_options('sample_by'), ['dittoFreqPlot']), false],
      plots: !is_ditto ? ['Data Representations', ['violin', 'box']] : ['Data Representations', ['vlnplot', 'boxplot', 'jitter', 'ridgeplot']],
      vars_use: ['Primary Data values to display (blank = all)', full_data, var_, 'Primary Data', discrete],
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
      size: !is_ditto() ? ['Point Size', 0.1, 50, undefined] : ['Point Size', 0.1, 25, 0.1],
      scale_by: [
        'Scale Y by counts or fraction',
        ['counts', 'fraction'],
        true,
        200,
        !is_ditto() ? x_by == y_by : false
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
      group_order: ['Reorder X-Axis Groupings?', full_data, group_by, discrete],
      var_order: ['Reorder Y-Axis / Y-Value Groupings?', full_data, var_, discrete],
      reduction_setup: [
        ['Dimensionality Reduction (DR)', 'x-axis DR Compenent #', 'y-axis DR Component #'],
        reduction_opts
      ],
      do_hover: ['Output as interactive (not flat) image?'],
      vlnplot_lineweight: ['Lineweight', 0, 1.5, 0.05],
      vlnplot_width: ['Width', 0.025, 1, 0.025],
      vlnplot_scaling: ['Volume Scaling Method', ['area', 'width', 'count'], false],
      boxplot_width: ['Width', 0.025, 0.5, 0.025],
      boxplot_color: ['Line Color'],
      boxplot_fill: ['Fill with "Color Data" colors?'],
      boxplot_lineweight: ['Lineweight', 0, 1.5, 0.05],
      jitter_size: ['Point Size', 0.05, 1, 0.05],
      jitter_width: ['Width', 0.025, 0.5, 0.025],
      jitter_color: ['Point Color'],
      ridgeplot_lineweight: ['Lineweight', 0, 1.5, 0.05],
      do_label: ['On/Off (discrete color data & flat output style only)'],
      labels_highlight: ['Highlight labels with a box?'],
      labels_repel: ['Allow movement to reduce overlaps?'],
      do_contour: ['Add density-based contours?'],
      do_ellipse: ['Surround groups in median-centered ellipses?'],
      legend_show: ['Include the legend?']
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
    group_by,
    var_,
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
  plots: MultiselectPiece,
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
  plot_subtitle: stringPiece,
  legend_title: stringPiece,
  xlab: stringPiece,
  ylab: stringPiece,
  x_by: nestedDropdownPiece,
  y_by: nestedDropdownPiece,
  var: nestedDropdownPiece,
  group_by: nestedDropdownPiece,
  color_by: nestedDropdownPiece,
  sample_by: nestedDropdownPiece,
  plots: MultiselectPiece,
  vars_use: MultiselectAfterDataChoicePiece_forDitto,
  color_order: ReorderPiece,
  size: sliderPiece,
  scale_by: dropdownPiece,
  x_scale: dropdownPiece,
  y_scale: dropdownPiece,
  cells_use: subsetDataFramePiece,
  group_order: ReorderCustomOnlyPiece,
  var_order: ReorderCustomOnlyPiece,
  reduction_setup: reductionSetupPiece,
  do_hover: checkboxPiece,
  vlnplot_lineweight: sliderPiece,
  vlnplot_width: sliderPiece,
  vlnplot_scaling: dropdownPiece,
  boxplot_width: sliderPiece,
  boxplot_color: stringPiece,
  boxplot_fill: checkboxPiece,
  boxplot_lineweight: sliderPiece,
  jitter_size: sliderPiece,
  jitter_width: sliderPiece,
  jitter_color: stringPiece,
  ridgeplot_lineweight: sliderPiece,
  do_label: checkboxPiece,
  labels_highlight: checkboxPiece,
  labels_repel: checkboxPiece,
  do_contour: checkboxPiece,
  do_ellipse: checkboxPiece,
  legend_show: checkboxPiece
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
  redefaults: DataEnvelope<DataEnvelope<DataEnvelope<any>>>,
  input_sets: DataEnvelope<DataEnvelope<string[]>>,
  components: DataEnvelope<Function>,
  plot_relabels?: DataEnvelope<string> | undefined
) {
  const preset = useMemo(() => data && data['preset'], [data]);
  const hide = useMemo(() => preset && Object.keys(preset), [preset]);
  const defaultValue = whichDefaults(setPlotType, preset, defaults, redefaults, input_sets);
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
  const discrete_columns: string[] = useMemo(() => {
    if (data == null) return [];
    if (data['discrete_cols'] == null) return df_columns; // Should build a warning here instead?
    return data['discrete_cols'];
  }, [data]);

  const columns: string[] | nestedOptionSet = (data != null && data['all_cols'] != null) ? data['all_cols'] : df_columns

  const reduction_opts: DataEnvelope<number[]> | null = useMemo(() => {
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
  const group_by =
    value && Object.keys(value).includes('group_by')
      ? (value.group_by as string | null)
      : null;
  const var_ =
    value && Object.keys(value).includes('var')
      ? (value.var as string | null)
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
    group_by,
    var_,
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
    onChange(some(whichDefaults(newType, preset, defaults, redefaults, input_sets)));
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
          plot_relabels == undefined ? Object.keys(input_sets) : plot_relabels,
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

  console.log(props.value);

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
  return VisualizationUI({data, onChange, value}, 'scatter_plot', defaults_plotly, redefaults_plotly, input_sets_plotly, components_plotly);
}

export function BarPlotly({
  data,
  onChange,
  value
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  return VisualizationUI({data, onChange, value}, 'bar_plot', defaults_plotly, redefaults_plotly, input_sets_plotly, components_plotly);
}

export function YPlotly({
  data,
  onChange,
  value
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  return VisualizationUI({data, onChange, value}, 'y_plot', defaults_plotly, redefaults_plotly, input_sets_plotly, components_plotly);
}

export function AnyPlotly({
  data,
  onChange,
  value
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  return VisualizationUI({data, onChange, value}, null, defaults_plotly, redefaults_plotly, input_sets_plotly, components_plotly, plot_relabels_plotly);
}

export function DittoDimPlot({
  data,
  onChange,
  value
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  return VisualizationUI({data, onChange, value}, 'dittoDimPlot', defaults_dittoseq, redefaults_dittoseq, input_sets_dittoseq, components_dittoseq);
}

export function DittoScatterPlot({
  data,
  onChange,
  value
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  return VisualizationUI({data, onChange, value}, 'dittoScatterPlot', defaults_dittoseq, redefaults_dittoseq, input_sets_dittoseq, components_dittoseq);
}

export function DittoBarPlot({
  data,
  onChange,
  value
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  return VisualizationUI({data, onChange, value}, 'dittoBarPlot', defaults_dittoseq, redefaults_dittoseq, input_sets_dittoseq, components_dittoseq);
}

export function DittoPlot({
  data,
  onChange,
  value
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  return VisualizationUI({data, onChange, value}, 'dittoPlot', defaults_dittoseq, redefaults_dittoseq, input_sets_dittoseq, components_dittoseq);
}

export function DittoFreqPlot({
  data,
  onChange,
  value
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  return VisualizationUI({data, onChange, value}, 'dittoFreqPlot', defaults_dittoseq, redefaults_dittoseq, input_sets_dittoseq, components_dittoseq);
}

export function AnyDittoSeq({
  data,
  onChange,
  value
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  return VisualizationUI({data, onChange, value}, null, defaults_dittoseq, redefaults_dittoseq, input_sets_dittoseq, components_dittoseq, plot_relabels_dittoseq);
}
