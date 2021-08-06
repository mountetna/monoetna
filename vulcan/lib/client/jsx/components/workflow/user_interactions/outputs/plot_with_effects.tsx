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
  const [updatedData, setUpdatedData] = useState([] as any);
  const [highlightedTrace, setHighlightedTrace] = useState(-1);
  const plotRef = useRef<HTMLDivElement>(null);

  const emphasizeTrace = (trace: any): any => {
    return {...trace};
  };

  const deemphasizeTrace = (trace: any): any => {
    return {
      ...trace,
      marker: {
        ...trace.marker,
        opacity: 0.1
      }
    };
  };

  useEffect(() => {
    if (-1 !== highlightedTrace) {
      let copy = [...data];
      for (var i = 0; i < copy.length; i++) {
        if (i === highlightedTrace) {
          copy[i] = emphasizeTrace(copy[i]);
        } else {
          copy[i] = deemphasizeTrace(copy[i]);
        }
      }
      setUpdatedData(copy);
    } else {
      setUpdatedData(data);
    }
  }, [data, highlightedTrace]);

  useLegendHover(
    plotRef.current,
    setHighlightedTrace
  );

  return (
    <div ref={plotRef}>
      <Plot
        data={updatedData}
        layout={layout}
        config={config}
      />
    </div>
  );
};

export default PlotWithEffects;
