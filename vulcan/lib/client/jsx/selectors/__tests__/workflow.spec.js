import XYPlot from 'etna-js/plots/components/xy_plot/xy_plot';

import {plotDataForStep} from '../workflow_selector';

describe('Workflow Selector', () => {
  it('correctly returns plot and data for XY plot', () => {
    const step = {
      name: 'ui_plot',
      run: 'xy.cwl',
      in: {
        series0: 'umap/umap_data',
        group0: 'umap/expression_matrix',
        series0__type: 'scatter',
        series1: 'umap/pca_data',
        series1__type: 'scatter'
      },
      out: []
    };

    const consignment = {
      umap_data: {
        rows: [
          [1, 2],
          [2, 1]
        ],
        num_rows: 2,
        num_cols: 2,
        col_names: ['x', 'y'],
        row_names: ['test1', 'test2']
      },
      pca_data: {
        rows: [
          [-10, 22],
          [24, 11]
        ],
        num_rows: 2,
        num_cols: 2,
        col_names: ['x', 'y'],
        row_names: ['pca1', 'pca2']
      },
      expression_matrix: {
        rows: [
          [21.3, -40.2, 99.9],
          [-100.1, 42.2, 74.21]
        ],
        num_rows: 2,
        num_cols: 3,
        col_names: ['gene1', 'gene2', 'gene3'],
        row_names: ['test1', 'test2']
      }
    };

    const parentWidth = 40;

    let {plot, data} = plotDataForStep(step, consignment, parentWidth);

    expect(plot.component).toEqual(XYPlot);
    expect(plot.plot_type).toEqual('xy');

    let {
      configuration: {layout}
    } = plot;

    expect(layout.width).toEqual(parentWidth);
    expect(layout.height).toEqual(30);

    expect(data.plot_series.length).toEqual(2);

    expect(data.plot_series.map((s) => s.name)).toEqual([
      'umap_data',
      'pca_data'
    ]);
    expect(data.plot_series.map((s) => s.series_type)).toEqual([
      'scatter',
      'scatter'
    ]);

    expect(data.plot_series[0].variables.x.values).toEqual([1, 2]);
    expect(data.plot_series[0].variables.y.values).toEqual([2, 1]);
    expect(data.plot_series[0].variables.label.values).toEqual([
      'test1',
      'test2'
    ]);
    expect(data.plot_series[1].variables.x.values).toEqual([-10, 24]);
    expect(data.plot_series[1].variables.y.values).toEqual([22, 11]);
    expect(data.plot_series[1].variables.label.values).toEqual([
      'pca1',
      'pca2'
    ]);
  });
});
