import React, {useContext, useRef, useEffect, useState} from 'react';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {VulcanContext} from '../../../contexts/vulcan';

import {selectConsignment} from 'etna-js/plots/selectors/consignment_selector';
import XYPlot from 'etna-js/plots/components/xy_plot/xy_plot';
import CategoryPlot from 'etna-js/plots/components/category_plot/category_plot';

import {plotDataForStep} from '../../../selectors/workflow_selector';
import {plotData} from '../../../../../../../etna/packages/etna-js/plots/plot_script';

export default function Plot({md5sum}) {
  const parentRef = useRef();
  const [parentDims, setParentDims] = useState({width: 0, height: 0});
  const [plot, setPlot] = useState(null);
  const [data, setData] = useState(null);

  useEffect(() => {
    console.log(parentRef);
    if (parentRef.current) {
      setParentDims({
        width: parentRef.current.offsetWidth,
        height: parentRef.current.offsetHeight
      });
    }
  }, [parentRef]);

  const {workflow, pathIndex, stepIndex, setData, setWorkflow} = useContext(
    VulcanContext
  );

  // The md5sum may come from the workflow step at some point?
  let browserState = useReduxState(
    browserStateOf(md5sum, workflow.steps[0][5])
  );

  useEffect(() => {
    browserState = useReduxState(
      browserStateOf(md5sum, workflow.steps[0][5], parentDims)
    );
  }, [parentDims]);

  let {consignment} = browserState;

  console.log(plot);
  console.log(data);
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

function browserStateOf(md5sum, step) {
  return (state) => {
    console.log(step);

    let consignment = selectConsignment(state, md5sum);

    if (!consignment) return null;

    return {consignment};
  };
}
