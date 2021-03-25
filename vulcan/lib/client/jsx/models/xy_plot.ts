import XYPlot from 'etna-js/plots/components/xy_plot/xy_plot';
import Plot from './plot';
import Domain from './domain';
import Series from './series';
import {WorkflowStep} from "../api_types";

// export default class XYPlotModel extends Plot {
//   public type: string = 'xy';
//   public component = XYPlot;
//
//   constructor(step: WorkflowStep, consignment: any, parentWidth: number) {
//     super(step, consignment, parentWidth);
//     this.type = 'xy';
//     this.component = XYPlot;
//
//     this.configObj = {
//       plot_type: this.type,
//       configuration: {
//         layout: {},
//         variables: {
//           x_label: 'x', // ??
//           y_label: 'y', // ??
//           x_min: null,
//           y_min: null,
//           x_max: null,
//           y_max: null
//         },
//         xLabel: 'x', // ??
//         yLabel: 'y', // ??
//         plot_series: [], // I think this is only used to pull data out of the consignment
//         access: 'view',
//         configuration: {
//           config: {}
//         }
//       }
//     };
//     this.dataObj = {
//       xdomain: new Domain([0]).asArray,
//       ydomain: new Domain([0]).asArray,
//       plot_series: []
//     };
//
//     this.defineSeries();
//     this.calculatePlotObj();
//     this.calculateLayout();
//     this.calculateDataObj();
//   }
//
//   defineSeries() {
//     let seriesRegex = /^series(?<index>\d+)$/;
//     let groupRegex = /^group(?<index>\d+)$/;
//     let seriesTypeRegex = /^series(?<index>\d+)__type$/;
//
//     // Maximum number of series is # of inputs / 2,
//     //   since each step input needs at least series# and
//     //   series#__type.
//     let inputNames = Object.keys(this.step.in);
//     this.series = new Array(Math.round(inputNames.length / 2)).fill(null);
//
//     // First, let's get all the series names
//     //   and types set up.
//     inputNames.forEach((input) => {
//       let indexMatch;
//       let name;
//       let type;
//       let xValues;
//       let yValues;
//       let labels;
//       let secondaryData;
//       let inputMatch;
//
//       if ((inputMatch = input.match(seriesRegex)) && inputMatch.groups) {
//         indexMatch = inputMatch.groups.index;
//         name = this.step.in[input].split('/')[1];
//
//         let consignmentData = this.consignment[name];
//
//         if (consignmentData) {
//           xValues = consignmentData.rows.map((r: any) => r[0]);
//           yValues = consignmentData.rows.map((r: any) => r[1]);
//           labels = consignmentData.row_names;
//         }
//       } else if ((inputMatch = input.match(groupRegex)) && inputMatch.groups) {
//         indexMatch = inputMatch.groups.index;
//         const secondaryDataName = this.step.in[input].split('/')[1];
//
//         secondaryData = this.consignment[secondaryDataName];
//       } else if ((inputMatch = input.match(seriesTypeRegex)) && inputMatch.groups) {
//         indexMatch = inputMatch.groups.index;
//         type = this.step.in[input];
//       }
//
//       if (!indexMatch) return;
//
//       const indexValue = parseInt(indexMatch, 10);
//
//       if (!this.series[indexValue]) {
//         this.series[indexValue] = new Series();
//       }
//
//       if (name) this.series[indexValue].setName(name);
//       if (type) this.series[indexValue].setType(type);
//       if (xValues) this.series[indexValue].setXValues(xValues);
//       if (yValues) this.series[indexValue].setYValues(yValues);
//       if (labels) this.series[indexValue].setLabels(labels);
//       if (secondaryData)
//         this.series[indexValue].setSecondaryData(secondaryData);
//     });
//   }
//
//   calculatePlotObj() {
//     this.configObj.configuration.variables = {
//       ...this.configObj.configuration.variables,
//       ...this.minsMaxes
//     };
//   }
//
//   calculateDataObj() {
//     let minMax = {...this.minsMaxes};
//     let xDomain = new Domain([minMax.x_min, minMax.x_max]);
//     let yDomain = new Domain([minMax.y_min, minMax.y_max]);
//
//     this.dataObj.xdomain = xDomain.asArray;
//     this.dataObj.ydomain = yDomain.asArray;
//
//     this.updateDataPlotSeries();
//   }
//
//   updateDataPlotSeries() {
//     this.dataObj.plot_series = this.validSeries.map((s) => s.asObj);
//   }
//
//   calculateLayout() {
//     // Calculate the plot layout given the dimensions
//     //   of a parent container.
//     this.configObj.configuration.layout = {
//       height: Math.round((3 * this.parentWidth) / 4),
//       width: this.parentWidth,
//       margin: {
//         bottom: 100,
//         left: 100,
//         right: 100,
//         top: 100
//       }
//     };
//   }
//
//   colorSeriesBy(seriesNumber: number, colorBy: string) {
//     if (!this.series[seriesNumber]) return;
//
//     this.series[seriesNumber].setColorBy(colorBy);
//
//     this.updateDataPlotSeries();
//   }
//
//   getSeriesColorOptions(seriesNumber: number) {
//     if (!this.series[seriesNumber] || !this.series[seriesNumber].secondaryData)
//       return null;
//
//     const {secondaryData} = this.series[seriesNumber];
//     return secondaryData ? secondaryData.row_names : [];
//   }
// }
