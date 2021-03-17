import React from 'react';
import Plotly from 'react-plotly.js';
import {autoColors} from 'etna-js/utils/colors';

const Plot = ({data, layout}) => {
  if (!data ||!layout) return null;

  console.log({layout, data});
  
  return <Plotly 
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
        width: 1024,
        height: 1024,
        paper_bgcolor: '#fff',
        plot_bgcolor: '#f0f0f0'
      }
    }
    config={
      {
        displaylogo: false,
        toImageButtonOptions: {
          format: 'svg',
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
};

export default Plot;
