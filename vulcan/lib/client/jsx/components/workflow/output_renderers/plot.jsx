import {useReduxState} from 'etna-js/hooks/useReduxState';

import {plotData} from 'etna-js/plots/plot_script';

export default function Plot({md5sum}) {
  const browserState = useReduxState(browserStateOf({md5sum}));

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

function browserStateOf({md5sum}) {
  return (state) => {
    // Where can we get `plot` from??
    let plot = {
      plot_type: 'xy',
      configuration: {
        layout: {},
        variables: {}
      }
    };
    let consignment = selectConsignment(state, md5sum);
    let data = consignment ? plotData(plot, consignment) : null;
    return {plot, data};
  };
}
