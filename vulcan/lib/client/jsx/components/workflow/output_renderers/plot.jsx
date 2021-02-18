import React, {useContext, useEffect, useState} from 'react';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {VulcanContext} from '../../../contexts/vulcan';

import Dropdown from 'etna-js/components/inputs/dropdown';
import {selectConsignment} from 'etna-js/plots/selectors/consignment_selector';

import CategoryPlot from 'etna-js/plots/components/category_plot/category_plot';

import {plotModelForStep} from '../../../selectors/workflow_selector';
import PlotColorPicker from './plot_color_picker';

export default function Plot({md5sum, parentRef}) {
  const [plot, setPlot] = useState(null);

  const {workflow, pathIndex, stepIndex} = useContext(VulcanContext);

  // The md5sum may come from the workflow step at some point?
  let browserState = useReduxState(browserStateOf(md5sum));

  let {consignment} = browserState;

  useEffect(() => {
    if (parentRef.current) {
      let parentWidth = parentRef.current.offsetWidth;
      let plotModel = plotModelForStep(
        workflow.steps[0][5], // use pathIndex and stepIndex once we have those tied in
        consignment,
        parentWidth
      );
      setPlot(plotModel);
    }
  }, [parentRef.current]);

  if (!plot) return null;

  let {
    configuration: {layout, variables}
  } = plot.config;
  let PlotComponent = plot.component;

  if (!PlotComponent) return null;

  return (
    <div>
      {plot.hasColorableSeries ? (
        <PlotColorPicker plot={plot}></PlotColorPicker>
      ) : (
        <PlotComponent
          data={plot.data}
          plot={plot.config}
          config_variables={variables}
          layout={layout}
          parent_width={plot.parentWidth}
        />
      )}
    </div>
  );
}

function browserStateOf(md5sum) {
  return (state) => {
    let consignment = selectConsignment(state, md5sum);

    return {consignment};
  };
}
