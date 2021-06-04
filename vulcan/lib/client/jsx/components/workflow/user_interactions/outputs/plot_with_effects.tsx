import React, {useRef} from 'react';
import Plot from 'react-plotly.js';
import {useLegendHover} from './plot_effects';

const PlotWithEffects = ({
  data,
  layout,
  config
}: {
  data: any;
  layout: any;
  config: any;
}) => {
  const plotRef = useRef<HTMLDivElement>(null);

  useLegendHover(plotRef);

  return (
    <div ref={plotRef}>
      <Plot data={data} layout={layout} config={config} />
    </div>
  );
};

export default PlotWithEffects;
