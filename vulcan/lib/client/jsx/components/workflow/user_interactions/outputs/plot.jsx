import React from 'react';
import Plot from 'react-plotly.js';
import {autoColors} from 'etna-js/utils/colors';

const PlotOutput = ({data: plots}) => {
  return <React.Fragment>
    {Object.keys(plots).map(k => {
        const plot = plots[k];
        if (!plot) return null;
        let { data, layout } = plot;
        if (!data ||!layout) return null;

        return <React.Fragment key={k}>
          <Plot
            data={data}
            layout={
              {
                ...layout,
                colorway: autoColors(40),
                xaxis: {
                  ...layout.xaxis,
                  gridcolor: '#ccc',
                  gridwidth: 1,
                  zerolinecolor: '#ccc',
                  zerolinewidth: 1,
                  linecolor: '#000'
                },
                yaxis: {
                  ...layout.yaxis,
                  gridcolor: '#ccc',
                  gridwidth: 1,
                  zerolinecolor: '#ccc',
                  zerolinewidth: 1,
                  linecolor: '#000'
                },
                width: 800,
                height: 800,
                paper_bgcolor: '#fff',
                plot_bgcolor: '#f0f0f0'
              }
            }
            config={
              {
                displaylogo: false,
                toImageButtonOptions: {
                  format: 'png',
                  filename: 'plot',
                  height: 500,
                  width: 700,
                  scale: 1
                },
                modeBarButtonsToRemove: [
                  'pan2d',
                  'select2d',
                  'zoomIn2d',
                  'zoomOut2d',
                  'autoScale2d',
                  'hoverClosestCartesian',
                  'hoverCompareCartesian',
                  'toggleSpikelines'
                ]
              }
            }
          />
        </React.Fragment>
      })}
  </React.Fragment>
};

export default PlotOutput;
