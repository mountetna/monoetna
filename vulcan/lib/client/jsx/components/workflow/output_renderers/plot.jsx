import React, {useContext} from 'react';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {VulcanContext} from '../../../contexts/vulcan';

import {plotData} from 'etna-js/plots/plot_script';
import {selectConsignment} from 'etna-js/plots/selectors/consignment_selector';
import XYPlot from 'etna-js/plots/components/xy_plot/xy_plot';
import CategoryPlot from 'etna-js/plots/components/category_plot/category_plot';
import Vector from 'etna-js/plots/models/vector';

export default function Plot({md5sum}) {
  const {workflow, pathIndex, stepIndex, setData, setWorkflow} = useContext(
    VulcanContext
  );

  // Will this be consistent?
  const primaryOutputName = workflow.steps[0][4].out[0];

  const browserState = useReduxState(browserStateOf(md5sum, primaryOutputName));

  let {plot, data, width} = browserState;

  if (!plot || !data) return null;

  let {
    plot_type,
    configuration: {layout, variables}
  } = plot;
  let PlotComponent;
  switch (plot_type) {
    case 'xy':
      PlotComponent = XYPlot;
      break;
    case 'category':
      PlotComponent = CategoryPlot;
      break;
    default:
      return null;
  }

  return (
    <PlotComponent
      data={data}
      plot={plot}
      config_variables={variables}
      layout={layout}
      parent_width={800}
    />
  );
}

function padInt(arrayValues, func) {
  const margin = 1.2;
  return Math.round(margin * func(...arrayValues));
}

function vectorize(values) {
  return new Vector(
    values.map((val) => ({
      label: null,
      value: val
    }))
  );
}

function browserStateOf(md5sum, primaryOutputName) {
  return (state) => {
    let consignment = selectConsignment(state, md5sum);

    if (!consignment) return {};

    let primaryConsignment = consignment[primaryOutputName];

    if (!primaryConsignment) return {};
    // Where can we get `plot` from?? -- possible to infer from
    //   the consignment?
    const x = primaryConsignment.col_names[0];
    const y = primaryConsignment.col_names[1];
    let xValues = primaryConsignment.rows.map((r) => r[0]);
    let yValues = primaryConsignment.rows.map((r) => r[1]);
    let plot = {
      plot_type: 'xy',
      configuration: {
        layout: {
          height: 400,
          width: 600,
          margin: {
            bottom: 100,
            left: 100,
            top: 50,
            right: 50
          }
        },
        variables: {
          x_label: x,
          y_label: y,
          x_min: Math.min(...xValues),
          y_min: Math.min(...yValues),
          x_max: Math.max(...xValues),
          y_max: Math.max(...yValues)
        },
        xLabel: x,
        yLabel: y,
        plot_series: [
          {
            name: primaryOutputName,
            series_type: 'scatter',
            variables: {
              x: x,
              y: y
            }
          }
        ],
        access: 'view',
        configuration: {
          config: {}
        }
      }
    };

    // let mockConsignment = {
    //   x: {
    //     vector: xValues.map((xVal) => ({
    //       label: null,
    //       value: xVal
    //     }))
    //   },
    //   y: {
    //     vector: yValues.map((yVal) => ({
    //       label: null,
    //       value: yVal
    //     }))
    //   },
    //   series0___x: {
    //     vector: xValues.map((xVal) => ({
    //       label: null,
    //       value: xVal
    //     }))
    //   },
    //   series0___y: {
    //     vector: yValues.map((yVal) => ({
    //       label: null,
    //       value: yVal
    //     }))
    //   },
    //   xy____xdomain: {
    //     vector: [padInt(xValues, Math.min), padInt(xValues, Math.max)]
    //   },
    //   xy____ydomain: {
    //     vector: [padInt(yValues, Math.min), padInt(yValues, Math.max)]
    //   }
    // };

    // let test = plotData(plot, mockConsignment);
    let labels = primaryConsignment.row_names;

    console.log(xValues);
    console.log(yValues);
    console.log(labels);

    let data = {
      xdomain: [padInt(xValues, Math.min), padInt(xValues, Math.max)],
      ydomain: [padInt(yValues, Math.min), padInt(yValues, Math.max)],
      plot_series: [
        {
          name: primaryOutputName,
          series_type: 'scatter',
          variables: {
            x: vectorize(xValues),
            y: vectorize(yValues),
            labels: vectorize(labels)
          }
        }
      ]
    };
    console.log(data);
    return {plot, data};
  };
}
