import XYPlot from 'etna-js/plots/components/xy_plot/xy_plot';
import Plot from './plot';
import Domain from './domain';
import Series from './series';

export default class XYPlotModel extends Plot {
  constructor(step, consignment, parentWidth) {
    super(step, consignment, parentWidth);
    this.type = 'xy';

    this.plotObj = {
      plot_type: this.type,
      component: XYPlot,
      configuration: {
        layout: {},
        variables: {
          x_label: 'x', // ??
          y_label: 'y', // ??
          x_min: null,
          y_min: null,
          x_max: null,
          y_max: null
        },
        xLabel: 'x', // ??
        yLabel: 'y', // ??
        plot_series: [], // I think this is only used to pull data out of the consignment
        access: 'view',
        configuration: {
          config: {}
        }
      },
      parentWidth: this.parentWidth
    };
    this.dataObj = {
      xdomain: new Domain([0]).asArray,
      ydomain: new Domain([0]).asArray,
      plot_series: []
    };

    this.defineSeries();
    this.calculatePlotObj();
    this.calculateLayout();
    this.calculateDataObj();
  }

  get validSeries() {
    // Because we may have null series in the List,
    //   we'll need to filter those out and only
    //   use the "real" series.
    return this.series.filter((s) => null !== s);
  }

  get minsMaxes() {
    // Calculate this across all series
    let current = {
      x_min: Math.min(), // Counterintuitively, returns Infinity
      x_max: Math.max(), // Counterintuitively, returns -Infinity
      y_min: Math.min(),
      y_max: Math.max()
    };
    this.validSeries.forEach((series) => {
      current.x_min = Math.min(...[...series.xValues, current.x_min]);
      current.x_max = Math.max(...[...series.xValues, current.x_max]);
      current.y_min = Math.min(...[...series.yValues, current.y_min]);
      current.y_max = Math.max(...[...series.yValues, current.y_max]);
    });
    return current;
  }

  defineSeries() {
    let seriesRegex = /^series(?<index>\d+)$/;
    let groupRegex = /^group(?<index>\d+)$/;
    let seriesTypeRegex = /^series(?<index>\d+)__type$/;

    // Maximum number of series is # of inputs / 2,
    //   since each step input needs at least series# and
    //   series#__type.
    let inputNames = Object.keys(this.step.in);
    this.series = new Array(Math.round(inputNames.length / 2)).fill(null);

    // First, let's get all the series names
    //   and types set up.
    inputNames.forEach((input) => {
      let indexMatch;
      let name;
      let type;
      let xValues;
      let yValues;
      let labels;
      if (input.match(seriesRegex)) {
        indexMatch = input.match(seriesRegex).groups.index;
        name = this.step.in[input].split('/')[1];

        let consignmentData = this.consignment[name];

        if (consignmentData) {
          xValues = consignmentData.rows.map((r) => r[0]);
          yValues = consignmentData.rows.map((r) => r[1]);
          labels = consignmentData.row_names;
        }
      } else if (input.match(groupRegex)) {
        indexMatch = input.match(groupRegex).groups.index;
        // Don't know how to add this in to the Series yet...
      } else if (input.match(seriesTypeRegex)) {
        indexMatch = input.match(seriesTypeRegex).groups.index;
        type = this.step.in[input];
      }

      if (!indexMatch) return;

      if (!this.series[indexMatch]) {
        this.series[indexMatch] = new Series();
      }

      if (name) this.series[indexMatch].setName(name);
      if (type) this.series[indexMatch].setType(type);
      if (xValues) this.series[indexMatch].setXValues(xValues);
      if (yValues) this.series[indexMatch].setYValues(yValues);
      if (labels) this.series[indexMatch].setLabels(labels);
    });
  }

  calculatePlotObj() {
    this.plotObj.configuration.variables = {
      ...this.plotObj.configuration.variables,
      ...this.minsMaxes
    };
  }

  calculateDataObj() {
    let minMax = {...this.minsMaxes};
    let xDomain = new Domain([minMax.x_min, minMax.x_max]);
    let yDomain = new Domain([minMax.y_min, minMax.y_max]);

    this.dataObj.xdomain = xDomain.asArray;
    this.dataObj.ydomain = yDomain.asArray;

    this.dataObj.plot_series = this.validSeries.map((s) => s.asObj);
  }

  calculateLayout() {
    // Calculate the plot layout given the dimensions
    //   of a parent container.
    this.plotObj.configuration.layout = {
      height: Math.round((3 * this.parentWidth) / 4),
      width: this.parentWidth,
      margin: {
        bottom: 100,
        left: 100,
        right: 100,
        top: 100
      }
    };
  }
}
