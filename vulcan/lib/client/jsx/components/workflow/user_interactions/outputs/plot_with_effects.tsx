import React, {useEffect, useRef, useState} from 'react';
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
  const [hoverData, setHoverData] = useState({} as any);
  const [legendHoverIndex, setLegendHoverIndex] = useState(-1);
  const plotRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    setHoverData(data);
  }, [data]);

  useLegendHover(plotRef.current, setLegendHoverIndex);

  return (
    <div ref={plotRef}>
      <Plot data={hoverData} layout={layout} config={config} />
    </div>
  );
};

export default PlotWithEffects;
