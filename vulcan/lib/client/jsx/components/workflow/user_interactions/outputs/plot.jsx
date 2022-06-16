import React from 'react';
import {autoColors} from 'etna-js/utils/colors';
import PlotWithEffects from './plot_with_effects';

export const PlotOutput = ({data}) => {
  console.log({data})
  return <React.Fragment>
    {Object.keys(data).map((k) => {
      if (!Object.keys(data[k]).includes('layout')) {
        return null
      } else {
        return PlotlyOutput(data)
      }
    })}
  </React.Fragment>
}

export const PlotlyOutput = ({data: plots}) => {
  return (
    <React.Fragment>
      {Object.keys(plots).map((k) => {
        const plot = plots[k];
        if (!plot) return null;
        let {data, layout} = plot;
        if (!data || !layout) return null;

        return (
          <PlotWithEffects
            key={k}
            data={data}
            layout={{
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
            }}
            config={{
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
                'lasso2d',
                'zoom2d',
                'autoScale2d',
                'hoverClosestCartesian',
                'hoverCompareCartesian',
                'toggleSpikelines'
              ]
            }}
          />
        );
      })}
    </React.Fragment>
  );
};
