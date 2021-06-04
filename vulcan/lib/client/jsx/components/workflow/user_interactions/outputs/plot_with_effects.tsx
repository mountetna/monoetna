import React, {useRef} from 'react';
import Plot from 'react-plotly.js';
import {useLegendHover} from './plot_effects';

const PlotWithEffects = ({data, layout}: {data: any; layout: any}) => {
  const plotRef = useRef<HTMLDivElement>(null);

  useLegendHover(plotRef);

  return (
    <div ref={plotRef}>
      <Plot data={data} layout={layout} />
    </div>
  );
};

export default PlotWithEffects;
