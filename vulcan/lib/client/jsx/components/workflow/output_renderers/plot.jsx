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
      parent_width={width}
    />
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
    let plot = {
      plot_type: 'xy',
      configuration: {
        layout: {
          height: 400,
          width: 600,
          margin: {}
        },
        variables: {
          x_label: x,
          y_label: y
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

    let xValues = primaryConsignment.rows.map((r) => r[0]);
    let yValues = primaryConsignment.rows.map((r) => r[1]);
    let labels = primaryConsignment.row_names;

    const margin = 1.2;
    console.log(xValues);
    console.log(yValues);

    let data = {
      xdomain: [margin * Math.min(xValues), margin * Math.max(xValues)],
      ydomain: [margin * Math.min(yValues), margin * Math.max(yValues)],
      plot_series: [
        {
          name: primaryOutputName,
          series_type: 'scatter',
          variables: {
            x: new Vector(xValues),
            y: new Vector(yValues),
            labels: new Vector(labels)
          }
        }
      ]
    };
    return {plot, data};
  };
}
