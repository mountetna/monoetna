import React from 'react';
import Plot from 'react-plotly.js';

export default function PlotlyOutput({data}) {
  if (!data || !data.data || !data.layout) return null;

  return <Plot data={data.data} layout={data.layout}></Plot>;
}
