import React, {Dispatch, PropsWithChildren, useMemo, useState} from 'react';
import * as _ from 'lodash';

import {DataEnvelope, OptionSet, WithInputParams} from '../../input_types';
import {useSetsDefault} from '../../useSetsDefault';
import {some} from '../../../../../selectors/maybe';
import {
  Accordion,
  AccordionDetails,
  AccordionSummary,
  Grid,
  Typography
} from '@material-ui/core';
import {pick} from 'lodash';
import {
  ReductionSetupPiece,
} from '../pieces/user_input_pieces';
import {ReorderCustomOnlyPiece, ReorderVizPiece} from '../pieces/reorder_piece';
import NestedDropdownMultiChoiceAdvancedPiece from '../pieces/nested_dropdown_multi_choice_advanced_piece';
import DropdownPiece, { DropdownPieceRct } from '../pieces/dropdown_piece';
import NestedDropdownPiece from '../pieces/nested_dropdown_piece';
import { OptionalSelectionDefinitionPiece } from '../pieces/grouping_pieces';
import { StringPiece } from '../pieces/string_piece';
import { CheckboxPiece } from '../pieces/checkbox_piece';
import { MultiselectStringAfterDataChoicePiece, MultiselectStringPiece } from '../pieces/multiselect_string_piece';
import { SliderPiece } from '../pieces/number_pieces';

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
    ?optional:
      continuous_cols:
        Needed in order to work most properly, but allowed to work without it. A vector of column names of data_frame which contain continuous data.
      discrete_cols:
        Needed in order to work most properly, but allowed to work without it. A vector of column names of data_frame which contain discrete data.
      all_cols:
        string[] or nestedOptionSet. Useful to give as nestedOptionSet when some organization to the potential column names is ideal, e.g. for metadata vs genes versus ADTs with single-cell data.
      reduction_opts:
        E.g. {<reduciton-name>: ['1', '2', ..., 'n']}. An object powering selection of dimensionality reduction, and then which components to use. Keys should be reduction names that exist in the single-cell object, and values should be all possible components to use for that dimesnaionality reduction. A special case exists where one key can be '_Recommended_'.
      preset:
        A hash where keys = input pieces that should be hidden from the user & values = the preset value that should be used for that input. E.g. The x_by and y_by inputs are hardset within the umap workflow as '0' and '1', so a user has no choice here and should not be able to adjust those fields!

  'value', object where (key,val) pairs mostly equate to (key,val) pairs that could be splatted into an archimedes visualization funciton.

  Major design notes:
    - Organized around having pairs of input-names and associated component-setups.
    - An 'input_sets' object defines the set of inputs that are used for a given visualization type.  Thus, after a plot-type choice, necessary inputs are known and the set of necessary ui pieces can be compiled.
    - (key,val) pairs of the optional 'presets' input will cause any inputs' component-setup to be removed from the displayed set, while also providing the value to give to 'value[key]'.

*/

function remove_hidden(
  vals: DataEnvelope<any>,
  hide: string[] | null | undefined
): DataEnvelope<any> {
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
    'primary features': ['reduction_setup', 'color_by', 'size', 'opacity'],
    titles: ['plot_title', 'plot_subtitle', 'legend_title', 'xlab', 'ylab'],
    'data focus': ['color_order', 'cells_use'],
    'addition: labels': ['do_label', 'labels_highlight', 'labels_repel'],
    'other additions': ['do_contour', 'do_ellipse'],
    'output style': ['do_hover', 'legend_show']
  },
  dittoScatterPlot: {
    'primary features': ['x_by', 'y_by', 'color_by', 'size', 'opacity'],
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
  dittoDotPlot: {
    'primary features': ['vars', 'group_by', 'scale'],
    titles: ['plot_title', 'plot_subtitle', 'legend_title', 'xlab', 'ylab'],
    'output style': ['legend_show'],
    'data focus': ['cells_use'],
  },
  dittoFreqPlot: {
    'primary features': ['var', 'vars_use', 'group_by', 'sample_by', 'scale_by', 'plots', 'color_by'],
    titles: ['plot_title', 'plot_subtitle', 'legend_title', 'xlab', 'ylab'],
    'vlnplot tweaks': ['vlnplot_lineweight', 'vlnplot_width', 'vlnplot_scaling'],
    'boxplot tweaks': ['boxplot_width', 'boxplot_color', 'boxplot_fill', 'boxplot_lineweight'],
    'jitter tweaks': ['jitter_size', 'jitter_width', 'jitter_color'],
    'ridgeplot tweaks': ['ridgeplot_lineweight'],
    'output style': ['legend_show', 'split_adjust_free_y'],
    'data focus': ['group_order', 'cells_use'],
  }
};

const plot_relabels_dittoseq: DataEnvelope<string> = {
  'dittoDimPlot: UMAP, PCA, etc.': 'dittoDimPlot',
  'dittoScatterPlot: e.g. gene x gene': 'dittoScatterPlot',
  'dittoBarPlot: compositional stacked bar plot': 'dittoBarPlot',
  'dittoPlot: violin, box, or ridge plot': 'dittoPlot',
  'dittoDotPlot: expression comparisons across cell groups': 'dittoDotPlot',
  'dittoFreqPlot: compositional violin, box, or ridge plot': 'dittoFreqPlot'
}

const defaults_dittoseq: DataEnvelope<any> = {
  x_by: null,
  y_by: null,
  var: null,
  vars: null,
  group_by: null,
  color_by: null,
  sample_by: 'make',
  plots: ['jitter', 'vlnplot'],
  scale: true,
  scale_by: 'fraction',
  size: 1,
  opacity: 1,
  plot_title: 'make',
  plot_subtitle: 'make',
  legend_title: 'make',
  xlab: 'make',
  ylab: 'make',
  color_order: 'unordered',
  group_order: 'make',
  var_order: 'make',
  vars_use: [],
  x_scale: 'linear',
  y_scale: 'linear',
  cells_use: false,
  reduction_setup: [null, '1', '2'],
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
  do_label: true,
  labels_highlight: true,
  labels_repel: true,
  do_contour: false,
  do_ellipse: false,
  legend_show: true,
  split_adjust_free_y: true
};

const redefaults_dittoseq: DataEnvelope<DataEnvelope<any>> = {
  dittoScatterPlot: {
    'color_by': 'make',
    'do_hover': false
  },
  dittoPlot: {
    'color_by': 'make'
  },
  dittoDimPlot: {
    'do_hover': false
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
    dittoDotPlot: {
      group_by: 'discrete',
      vars: 'continuous',
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
  options: OptionSet,
  full_data: DataEnvelope<any>,
  plot_type: string | null,
  continuous: OptionSet,
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
  function add_make(add_to: OptionSet, for_plot_types: string[]) {
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
      vars: ['Primary Data', get_options('vars')],
      group_by: ['Groupings Data (often the x-axis)', get_options('group_by'), false],
      color_by: [
        plot_type!=null && ['dittoPlot', 'dittoFreqPlot'].includes(plot_type) ?
          'Color Data (\'make\' = same as Groupings Data)' : 'Color Data',
        add_make(get_options('color_by'), ['scatter_plot', 'y_plot', 'dittoPlot', 'dittoFreqPlot', 'dittoScatterPlot']),
        false],
      sample_by: ['Sample Data (try \'_sc_seq_ids_\' if there; \'make\' = ignore)', add_make(get_options('sample_by'), ['dittoFreqPlot']), false],
      plots: !is_ditto ? [
        'Data Representations',
        ['box', 'violin'],
        ['box', 'violin'],
        ['box', 'violin']
      ] : [
        'Data Representations',
        ['vlnplot', 'boxplot', 'jitter', 'ridgeplot'],
        ['vlnplot', 'boxplot', 'jitter', 'ridgeplot'],
        ['jitter', 'vlnplot']
      ],
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
      opacity: ['Point Opacity', 0.05, 1, 0.05],
      scale_by: [
        'Scale Y by counts or fraction',
        ['counts', 'fraction'],
        true,
        200,
        !is_ditto() ? x_by == y_by : false
      ],
      scale: ['scale (z-score) summary values?'],
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
        options,
        false,
        'secondary'
      ],
      x_order: ['Order of X-Axis Groupings', full_data, x_by, discrete],
      y_order: ['Order of Y-Axis Groupings', full_data, y_by, discrete],
      group_order: ['Reorder X-Axis Groupings?', full_data, group_by, discrete],
      var_order: ['Reorder Y-Axis / Y-Value Groupings?', full_data, var_, discrete],
      reduction_setup: [
        ['Dimensionality Reduction (DR)', 'x-axis DR Compenent #', 'y-axis DR Component #'],
        reduction_opts
      ],
      do_hover: ['Output as interactive image? (automatically ignored if too many data points)'],
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
      legend_show: ['Include the legend?'],
      split_adjust_free_y: ['Allow y-axis to vary between facets?']
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
  plot_title: StringPiece,
  legend_title: StringPiece,
  xlab: StringPiece,
  ylab: StringPiece,
  x_by: DropdownPiece,
  y_by: DropdownPiece,
  color_by: DropdownPiece,
  plots: MultiselectStringPiece,
  color_order: ReorderVizPiece,
  order_when_continuous_color: CheckboxPiece,
  size: SliderPiece,
  scale_by: DropdownPiece,
  x_scale: DropdownPiece,
  y_scale: DropdownPiece,
  rows_use: OptionalSelectionDefinitionPiece,
  x_order: ReorderVizPiece,
  y_order: ReorderVizPiece
};

const components_dittoseq: DataEnvelope<Function> = {
  plot_title: StringPiece,
  plot_subtitle: StringPiece,
  legend_title: StringPiece,
  xlab: StringPiece,
  ylab: StringPiece,
  x_by: NestedDropdownPiece,
  y_by: NestedDropdownPiece,
  var: NestedDropdownPiece,
  vars: NestedDropdownMultiChoiceAdvancedPiece,
  group_by: NestedDropdownPiece,
  color_by: NestedDropdownPiece,
  sample_by: NestedDropdownPiece,
  plots: MultiselectStringPiece,
  vars_use: MultiselectStringAfterDataChoicePiece,
  color_order: ReorderVizPiece,
  size: SliderPiece,
  opacity: SliderPiece,
  scale_by: DropdownPiece,
  scale: CheckboxPiece,
  x_scale: DropdownPiece,
  y_scale: DropdownPiece,
  cells_use: OptionalSelectionDefinitionPiece,
  group_order: ReorderCustomOnlyPiece,
  var_order: ReorderCustomOnlyPiece,
  reduction_setup: ReductionSetupPiece,
  do_hover: CheckboxPiece,
  vlnplot_lineweight: SliderPiece,
  vlnplot_width: SliderPiece,
  vlnplot_scaling: DropdownPiece,
  boxplot_width: SliderPiece,
  boxplot_color: StringPiece,
  boxplot_fill: CheckboxPiece,
  boxplot_lineweight: SliderPiece,
  jitter_size: SliderPiece,
  jitter_width: SliderPiece,
  jitter_color: StringPiece,
  ridgeplot_lineweight: SliderPiece,
  do_label: CheckboxPiece,
  labels_highlight: CheckboxPiece,
  labels_repel: CheckboxPiece,
  do_contour: CheckboxPiece,
  do_ellipse: CheckboxPiece,
  legend_show: CheckboxPiece,
  split_adjust_free_y: CheckboxPiece
};

export const ComponentUse = ({
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
  const preset = useMemo(() => data && 'preset' in data ? data['preset'] : undefined, [data]);
  const hide = useMemo(() => preset && Object.keys(preset), [preset]);
  const defaultValue = whichDefaults(setPlotType, preset, defaults, redefaults, input_sets);
  const value = useSetsDefault(defaultValue, props.value, onChange, 'plot_setup');
  const [expandedDrawers, setExpandedDrawers] = useState(['primary features']);
  const plotType =
    value && value['plot_type'] ? (value['plot_type'] as string) : null;

  const data_frame: DataEnvelope<any> = useMemo(() => {
    return data['data_frame'];
  }, [data]);

  const df_columns: string[] = useMemo(() => {
    if (data_frame == null || 0 === Object.keys(data_frame).length) return [];
    return Object.keys(data_frame);
  }, [data_frame]);

  const continuous_columns: OptionSet = useMemo(() => { // NOTE: these data don't necessarily need to be contained within the given data_frame (to accomodate for visualizing genomics data)
    return 'continuous_opts' in data ? data['continuous_opts'] : df_columns;
  }, [data]);
  const discrete_columns: string[] = useMemo(() => {
    return 'discrete_opts' in data ? data['discrete_opts'] : df_columns;
  }, [data]);
  const columns: OptionSet = useMemo(() => {
    return 'all_opts' in data ? data['all_opts'] : df_columns;
  }, [data]);

  const reduction_opts: DataEnvelope<string[]> | null = useMemo(() => {
    return 'reduction_opts' in data ? data['reduction_opts'] : null;
  }, [data])

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

  const shownSetupValues: ReturnType<typeof remove_hidden> = useMemo(() => {
    if (plotType == null) return {};
    const initial = {...value};
    if (Object.keys(initial).includes('plot_type')) {
      delete initial['plot_type'];
    }
    return remove_hidden(initial, hide);
  }, [value, hide]);

  const updatePlotType = (newType: string | null, key?: string) => {
    onChange({plot_setup: some(whichDefaults(newType, preset, defaults, redefaults, input_sets))});
  };

  const updateValue = (newValue: any, key: string, prevValues = {...value}) => {
    prevValues[key] = newValue;
    onChange({plot_setup: some(prevValues)});
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
        {DropdownPiece(
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
          const group_values: DataEnvelope<any> = pick(shownSetupValues, val as string[]);
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
                  <div key={key} style={{paddingTop: '8px'}}>
                    <ComponentUse
                      key={key}
                      k={key}
                      value={val}
                      extra_inputs={extra_inputs[key]}
                      updateValue={updateValue}
                      comps={components}
                    />
                  </div>
                );
              })}
            </InputWrapper>
          ) : null;
        })}
      </Grid>
    );

  // console.log({
  //   VizValue: value,
  //   data: data
  // });

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

export function DittoDotPlot({
  data,
  onChange,
  value
}: WithInputParams<{}, DataEnvelope<any>, any>) {
  return VisualizationUI({data, onChange, value}, 'dittoDotPlot', defaults_dittoseq, redefaults_dittoseq, input_sets_dittoseq, components_dittoseq);
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
