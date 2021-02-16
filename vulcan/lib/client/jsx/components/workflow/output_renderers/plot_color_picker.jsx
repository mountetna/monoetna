import React, {useEffect, useState} from 'react';

import Dropdown from 'etna-js/components/inputs/dropdown';

export default function PlotColorPicker({plot}) {
  const [seriesIndex, setSeriesIndex] = useState(null);
  const [colorByIndex, setColorByIndex] = useState(null);
  const [showPlot, setShowPlot] = useState(null);

  useEffect(() => {
    if (null !== seriesIndex && null !== colorByIndex) {
      let colorBy = plot.getSeriesColorOptions(seriesIndex)[colorByIndex];
      plot.colorSeriesBy(seriesIndex, colorBy);
      setShowPlot(true);
    }
  }, [seriesIndex, colorByIndex]);

  function handleOnSelectSeries(index) {
    setSeriesIndex(index);
    setShowPlot(false);
  }

  function handleOnSelectColorBy(index) {
    setColorByIndex(index);
    setShowPlot(false);
  }

  let {
    configuration: {layout, variables}
  } = plot.config;
  let PlotComponent = plot.component;

  return (
    <div>
      <div>
        <p>Select a series to color by:</p>
        <Dropdown
          list={plot.validSeries.map((s) => s.name)}
          selected_index={seriesIndex}
          onSelect={handleOnSelectSeries}
        ></Dropdown>
      </div>
      {null !== seriesIndex ? (
        <div>
          <p>Select a data set to color by:</p>
          <Dropdown
            list={plot.getSeriesColorOptions(seriesIndex)}
            selected_index={colorByIndex}
            onSelect={handleOnSelectColorBy}
          ></Dropdown>
        </div>
      ) : null}
      {showPlot ? (
        <PlotComponent
          data={plot.data}
          plot={plot.config}
          config_variables={variables}
          layout={layout}
          parent_width={plot.parentWidth}
        />
      ) : null}
    </div>
  );
}
