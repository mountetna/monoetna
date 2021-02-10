import React, {useContext, useRef, useEffect, useState} from 'react';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {VulcanContext} from '../../../contexts/vulcan';

import {selectConsignment} from 'etna-js/plots/selectors/consignment_selector';
import XYPlot from 'etna-js/plots/components/xy_plot/xy_plot';
import CategoryPlot from 'etna-js/plots/components/category_plot/category_plot';
import Vector from 'etna-js/plots/models/vector';

export default function Plot({md5sum}) {
  const parentRef = useRef();
  const [parentDims, setParentDims] = useState({width: 1600, height: 1200});

  useEffect(() => {
    if (parentRef.current) {
      setParentDims({
        width: parentRef.current.offsetWidth,
        height: parentRef.current.offsetHeight
      });
    }
  }, []);

  const {workflow, pathIndex, stepIndex, setData, setWorkflow} = useContext(
    VulcanContext
  );

  // The md5sum may come from the workflow step at some point?
  const browserState = useReduxState(
    browserStateOf(md5sum, workflow.steps[0][5], parentDims)
  );

  let {plot, data} = browserState;

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
    <div ref={parentRef}>
      <PlotComponent
        data={data}
        plot={plot}
        config_variables={variables}
        layout={layout}
        parent_width={parentDims.width}
      />
    </div>
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

function browserStateOf(md5sum, step, dimensions) {
  return (state) => {
    // Will this be consistent?
    const primaryOutputName = step.out[0];

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

    let labels = primaryConsignment.row_names;

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
            label: vectorize(labels)
          }
        }
      ]
    };

    return {plot, data};
  };
}
