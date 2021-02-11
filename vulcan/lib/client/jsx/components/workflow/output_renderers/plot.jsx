import React, {useContext, useRef, useEffect, useState} from 'react';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {VulcanContext} from '../../../contexts/vulcan';

import {selectConsignment} from 'etna-js/plots/selectors/consignment_selector';

import CategoryPlot from 'etna-js/plots/components/category_plot/category_plot';

import {plotDataForStep} from '../../../selectors/workflow_selector';

export default function Plot({md5sum, parentRef}) {
  const [plotConfig, setPlotConfig] = useState(null);
  const [plotData, setPlotData] = useState(null);

  const {workflow, pathIndex, stepIndex} = useContext(VulcanContext);

  // The md5sum may come from the workflow step at some point?
  let browserState = useReduxState(browserStateOf(md5sum));

  let {consignment} = browserState;

  useEffect(() => {
    if (parentRef.current) {
      let parentWidth = parentRef.current.offsetWidth;
      let {plot, data} = plotDataForStep(
        workflow.steps[0][5], // use pathIndex and stepIndex once we have those tied in
        consignment,
        parentWidth
      );
      setPlotConfig(plot);
      setPlotData(data);
    }
  }, [parentRef.current]);

  if (!plotConfig || !plotData) return null;

  let {
    configuration: {layout, variables},
    parentWidth
  } = plotConfig;
  let PlotComponent = plotConfig.component;

  if (!PlotComponent) return null;

  return (
    <div ref={parentRef}>
      <PlotComponent
        data={plotData}
        plot={plotConfig}
        config_variables={variables}
        layout={layout}
        parent_width={parentWidth}
      />
    </div>
  );
}

function browserStateOf(md5sum) {
  return (state) => {
    let consignment = selectConsignment(state, md5sum);
    console.log(consignment);
    return {consignment};
  };
}
