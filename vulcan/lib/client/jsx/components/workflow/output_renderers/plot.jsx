import React, {useContext} from 'react';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {VulcanContext} from '../../../contexts/vulcan';

import {plotData} from 'etna-js/plots/plot_script';
import {selectConsignment} from 'etna-js/plots/selectors/consignment_selector';
import XYPlot from 'etna-js/plots/components/xy_plot/xy_plot';
import CategoryPlot from 'etna-js/plots/components/category_plot/category_plot';

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
          width: 600
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

    let data = {
      xdomain: [-500, 500],
      ydomain: [-500, 500],
      plot_series: [
        {
          name: primaryOutputName,
          series_type: 'scatter',
          variables: {
            x: primaryConsignment.rows.map((r) => r[0]),
            y: primaryConsignment.rows.map((r) => r[1])
          }
        }
      ]
    };
    return {plot, data};
  };
}
