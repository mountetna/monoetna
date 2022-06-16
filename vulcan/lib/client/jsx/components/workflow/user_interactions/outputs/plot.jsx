import React from 'react';
import {autoColors} from 'etna-js/utils/colors';
import PlotWithEffects from './plot_with_effects';

export const PlotOutput = ({data}) => {
  const plotlys = data.raw
  const pngs = Object.fromEntries(Object.entries(data.url).
    filter(([key, val]) => {
      return !Object.keys(data.raw).includes(key) && val.endsWith(".png")
    }))
  return <React.Fragment>
    {(Object.keys(plotlys).length>0) ? PlotlyOutput({data: plotlys}) : null}
    {(Object.keys(pngs).length>0) ? PngOutput({data: pngs}) : null}
  </React.Fragment>
}

export const PngOutput = ({data}) => {
  return <React.Fragment>
    {Object.values(data).map((url, ind) =>
        <React.Fragment key={url}>
          <img key={ind} src={url}/>
        </React.Fragment>
    )}
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
