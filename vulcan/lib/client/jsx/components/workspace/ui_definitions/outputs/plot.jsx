import React, {useContext, useState, useEffect} from 'react';
import {autoColors} from 'etna-js/utils/colors';
import PlotWithEffects from './plot_with_effects';
import {VulcanContext} from '../../../../contexts/vulcan_context';

export const PngOutput = ({data}) => {
  let {state, getImage} = useContext(VulcanContext);
  let {workspace, projectName} = state;
  const [urls, setUrls] = useState([]);
  const [pulling, setPulling] = useState(false);

  const pullImages = async (data) => {
    if (!workspace) {
      setPulling(false);
      return
    }
    const imgs = Object.keys(data)
    const newUrls = [];
    for (let ind in imgs) {
      const res = await getImage(projectName, workspace.workspace_id, imgs[ind]);
      const imageBlob = await res.blob();
      const url = URL.createObjectURL(imageBlob);
      newUrls.push(url);
    }
    setUrls(newUrls);
    setPulling(false);
  }

  useEffect(() => {
    if (!pulling && Object.values(data).length > urls.length) {
      setPulling(true);
      pullImages(data);
    }
  }, [data, urls])

  return <React.Fragment>
    {urls.length == 0 ? null : urls.map((val, ind) => {
      return <React.Fragment key={ind}>
          <img src={val} width="800"/>
        </React.Fragment>
    })}
  </React.Fragment>
}

export const PlotlyOutput = ({data: plots}) => {
  return (
    <React.Fragment>
      {Object.keys(plots).map((k) => {
        let plot = {...plots}[k];
        if (!plot) return null;
        let {data, layout} = plot;
        if (!data || !layout) return null; // Todo: make this a user visible warning telling them to reach out to us!

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

export const PlotOutput = ({data, url}) => {
  const plotlys = Object.fromEntries(
    Object.entries(data).filter(([key, val])=> {
      return typeof(val)==='object'
    })
  )
  const pngs = Object.fromEntries(
    Object.entries(data).filter(([key, val]) => {
      return typeof(val)==='string' && val.length > 0
    })
  )
  return <React.Fragment>
    {(Object.keys(plotlys).length>0) ? PlotlyOutput({data: plotlys}) : PngOutput({data: pngs})}
  </React.Fragment>
}
