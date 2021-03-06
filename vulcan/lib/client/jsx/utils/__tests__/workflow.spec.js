import Gradient from 'javascript-color-gradient';

import {autoColors} from 'etna-js/utils/colors';
import XYPlot from 'etna-js/plots/components/xy_plot/xy_plot';

import {
  plotModelForStep,
  validPath,
  validStep,
  defaultInputValues,
  allInputsDefined,
  uiStepOptions,
  missingUiInputs,
  inputNamesToHashStub
} from '..//workflow';

describe('Workflow Utils', () => {
  describe('XY Plots', () => {
    it('correctly returns plot model', () => {
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

      let plot = plotModelForStep(step, consignment, parentWidth);

      let plotConfig = plot.config;
      let plotData = plot.data;

      expect(plot.component).toEqual(XYPlot);
      expect(plot.type).toEqual('xy');

      let {
        configuration: {layout}
      } = plotConfig;

      expect(layout.width).toEqual(parentWidth);
      expect(layout.height).toEqual(30);

      expect(plotData.plot_series.length).toEqual(2);

      expect(plotData.plot_series.map((s) => s.name)).toEqual([
        'umap_data',
        'pca_data'
      ]);
      expect(plotData.plot_series.map((s) => s.series_type)).toEqual([
        'scatter',
        'scatter'
      ]);

      expect(plotData.plot_series[0].variables.x.values).toEqual([1, 2]);
      expect(plotData.plot_series[0].variables.y.values).toEqual([2, 1]);
      expect(plotData.plot_series[0].variables.label.values).toEqual([
        'test1',
        'test2'
      ]);
      expect(plotData.plot_series[1].variables.x.values).toEqual([-10, 24]);
      expect(plotData.plot_series[1].variables.y.values).toEqual([22, 11]);
      expect(plotData.plot_series[1].variables.label.values).toEqual([
        'pca1',
        'pca2'
      ]);
    });

    it('can color a series', () => {
      const step = {
        name: 'ui_plot',
        run: 'xy.cwl',
        in: {
          series0: 'umap/umap_data',
          group0: 'umap/expression_matrix',
          series0__type: 'scatter'
        },
        out: []
      };

      const consignment = {
        umap_data: {
          rows: [
            [1, 2],
            [2, 1],
            [3, 4]
          ],
          num_rows: 2,
          num_cols: 2,
          col_names: ['x', 'y'],
          row_names: ['test1', 'test2', 'test3']
        },
        expression_matrix: {
          rows: [
            [21.3, -100, 99.9],
            [-100.1, 100, 74.21],
            [5.12, 0, 51.89]
          ],
          num_rows: 2,
          num_cols: 3,
          col_names: ['test1', 'test2', 'test3'],
          row_names: ['gene1', 'gene2', 'gene3']
        }
      };

      const parentWidth = 40;

      let plot = plotModelForStep(step, consignment, parentWidth);
      expect(plot.hasColorableSeries).toEqual(true);

      let plotData = plot.data;

      expect(plotData.plot_series[0].variables.nodeColor).toEqual(null);

      const colorOptions = plot.getSeriesColorOptions(0);

      expect(colorOptions.length).toEqual(
        consignment.expression_matrix.col_names.length
      );

      plot.colorSeriesBy(0, colorOptions[1]);

      plotData = plot.data;

      let gene2Colors = plotData.plot_series[0].variables.nodeColor;

      expect(gene2Colors.length).toEqual(3);

      plot.colorSeriesBy(0, colorOptions[0]);

      plotData = plot.data;

      let gene1Colors = plotData.plot_series[0].variables.nodeColor;

      expect(gene1Colors.length).toEqual(3);
      expect(gene1Colors).not.toEqual(gene2Colors);
    });
  });

  describe('validPath', () => {
    it('returns false if invalid', () => {
      let results = validPath({
        workflow: {},
        pathIndex: 0
      });

      expect(results).toEqual(false);

      results = validPath({
        workflow: {
          steps: [[]]
        },
        pathIndex: 1
      });

      expect(results).toEqual(false);
    });

    it('returns true if invalid', () => {
      let results = validPath({
        workflow: {
          steps: [[]]
        },
        pathIndex: 0
      });

      expect(results).toEqual(true);
    });
  });

  describe('validStep', () => {
    it('returns false if invalid', () => {
      let results = validStep({
        workflow: {},
        pathIndex: 0,
        stepIndex: 0
      });

      expect(results).toEqual(false);

      results = validStep({
        workflow: {
          steps: [[]]
        },
        pathIndex: 0,
        stepIndex: 1
      });

      expect(results).toEqual(false);
    });

    it('returns true if invalid', () => {
      let results = validStep({
        workflow: {
          steps: [[1]]
        },
        pathIndex: 0,
        stepIndex: 0
      });

      expect(results).toEqual(true);
    });
  });

  describe('defaultInputValues', () => {
    it('returns values', () => {
      const workflow = {
        class: 'Workflow',
        cwlVersion: 'v1.1',
        inputs: {
          bool_input: {
            default: true,
            label: 'Sample boolean',
            type: 'boolean'
          },
          int_input: {
            default: 42,
            label: 'Sample input',
            type: 'int'
          },
          no_default: {
            default: null,
            type: 'int'
          }
        },
        outputs: {
          sample_data: {
            outputSource: 'final_step/sample_data',
            type: 'File'
          }
        },
        steps: []
      };

      let results = defaultInputValues(workflow);
      expect(results).toEqual({bool_input: true, int_input: 42});
    });
  });

  describe('allInputsDefined', () => {
    it('returns true when defined', () => {
      const workflow = {
        class: 'Workflow',
        cwlVersion: 'v1.1',
        inputs: {
          bool_input: {
            default: true,
            label: 'Sample boolean',
            type: 'boolean'
          },
          int_input: {
            default: 42,
            label: 'Sample input',
            type: 'int'
          },
          no_default: {
            default: null,
            type: 'int'
          }
        },
        outputs: {
          sample_data: {
            outputSource: 'final_step/sample_data',
            type: 'File'
          }
        },
        steps: []
      };

      expect(
        allInputsDefined(workflow, {
          int_input: 1,
          bool_input: true,
          no_default: 2
        })
      ).toEqual(true);
    });

    it('returns false if null or undefined', () => {
      const workflow = {
        class: 'Workflow',
        cwlVersion: 'v1.1',
        inputs: {
          bool_input: {
            default: true,
            label: 'Sample boolean',
            type: 'boolean'
          },
          int_input: {
            default: 42,
            label: 'Sample input',
            type: 'int'
          },
          no_default: {
            default: null,
            type: 'int'
          }
        },
        outputs: {
          sample_data: {
            outputSource: 'final_step/sample_data',
            type: 'File'
          }
        },
        steps: []
      };

      expect(
        allInputsDefined(workflow, {
          bool_input: true,
          int_input: null,
          no_default: 2
        })
      ).toEqual(false);

      expect(
        allInputsDefined(workflow, {
          bool_input: true,
          int_input: 1,
          no_default: undefined
        })
      ).toEqual(false);

      expect(
        allInputsDefined(workflow, {
          bool_input: true,
          int_input: 1,
          no_default: 3,
          user_input: null
        })
      ).toEqual(false);

      expect(
        allInputsDefined(workflow, {
          bool_input: true,
          int_input: 1,
          no_default: 3,
          user_input: NaN
        })
      ).toEqual(false);
    });

    it('returns false if key missing', () => {
      const workflow = {
        class: 'Workflow',
        cwlVersion: 'v1.1',
        inputs: {
          bool_input: {
            default: true,
            label: 'Sample boolean',
            type: 'boolean'
          },
          int_input: {
            default: 42,
            label: 'Sample input',
            type: 'int'
          },
          no_default: {
            default: null,
            type: 'int'
          }
        },
        outputs: {
          sample_data: {
            outputSource: 'final_step/sample_data',
            type: 'File'
          }
        },
        steps: []
      };

      expect(
        allInputsDefined(workflow, {
          bool_input: true,
          int_input: 1
        })
      ).toEqual(false);
    });
  });

  describe('uiStepOptions', () => {
    it('retrieves data from a required step as an Array', () => {
      let step = {
        name: 'foo',
        in: [
          {
            source: ['previous', 'output']
          }
        ]
      };

      const status = [
        [
          {
            name: 'previous',
            data: {
              output: 'blah'
            }
          },
          {
            name: 'alternate',
            data: {
              output: [1, 2, 3]
            }
          }
        ]
      ];

      let results = uiStepOptions({step, status, pathIndex: 0});
      expect(results).toEqual(['blah']);

      step = {
        name: 'foo',
        in: [
          {
            source: ['alternate', 'output']
          }
        ]
      };

      results = uiStepOptions({step, status, pathIndex: 0});
      expect(results).toEqual([1, 2, 3]);
    });

    it('returns empty list when step has no data', () => {
      const step = {
        name: 'foo',
        in: [
          {
            source: ['previous', 'output']
          }
        ]
      };

      let status = [
        [
          {
            name: 'previous'
          }
        ]
      ];

      let results = uiStepOptions({step, status, pathIndex: 0});
      expect(results).toEqual([]);

      status = [
        [
          {
            name: 'previous',
            data: {}
          }
        ]
      ];

      results = uiStepOptions({step, status, pathIndex: 0});
      expect(results).toEqual([]);
    });
  });

  describe('missingUiInputs', () => {
    it('returns input names not in the session inputs', () => {
      const step = {
        out: ['output'],
        name: 'step1'
      };

      const session = {
        inputs: {
          a: 123
        }
      };

      let results = missingUiInputs(step, session);
      expect(results).toEqual(['step1/output']);
    });

    it('does not return input names already in the session inputs', () => {
      const step = {
        out: ['output'],
        name: 'step1'
      };

      const session = {
        inputs: {
          'step1/output': 123
        }
      };

      let results = missingUiInputs(step, session);
      expect(results).toEqual([]);
    });
  });

  describe('inputNamesToHashStub', () => {
    it('reduces array of names to Hash with null values', () => {
      const names = ['step1/output1', 'step1/output2'];

      let result = inputNamesToHashStub(names);
      expect(result).toEqual({
        'step1/output1': null,
        'step1/output2': null
      });
    });
  });
});
