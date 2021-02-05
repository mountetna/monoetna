import React, {useState, useContext, useEffect} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from 'etna-js/actions/message_actions';

import {VulcanContext} from '../../contexts/vulcan';
import {getData} from '../../api/vulcan';

import Plot from './output_renderers/plot';
import Raw from './output_renderers/raw';
import {stepDataUrls, stepData} from '../../selectors/workflow_selector';

const loadingDiv = (
  <div className='browser'>
    <div id='loader-container'>
      <div className='loader'>Loading...</div>
    </div>
  </div>
);

export default function Output() {
  const invoke = useActionInvoker();
  const {workflow, pathIndex, stepIndex, setData} = useContext(VulcanContext);

  const [outputData, setOutputData] = useState(null);
  const [outputIndex, setOutputIndex] = useState(0);
  const [outputDatum, setOutputDatum] = useState(null);

  // useEffect(() => {
  //   let dataUrls = stepDataUrls(workflow, pathIndex);
  //   dataUrls.forEach((url) => {
  //     getData(url)
  //       .then((data) => {
  //         setData(url, data);
  //       })
  //       .catch((e) => {
  //         invoke(showMessages(e));
  //       });
  //   });
  // }, [pathIndex, stepIndex]);

  useEffect(() => {
    // When the workflow updates, let's pull out
    //   the data for the current step.
    setOutputData(stepData(workflow, pathIndex, stepIndex));
  }, [workflow]);

  useEffect(() => {
    setData('https://vulcan.test/blob', {
      x: [1, 2, 3, 4, 5],
      y: [5, 4, 3, 2, 1]
    });
  }, []);

  useEffect(() => {
    setOutputDatum(outputData[outputIndex]);
  }, [outputData, outputIndex]);

  const viewOptions = {
    plot: Plot,
    raw: Raw
  };

  let Component;
  function handleOnSelectView(index) {
    Component = viewOptions[Object.keys(viewOptions)[index]];
  }

  function handleOnSelectOutput(index) {
    setOutputIndex(index);
  }

  if (!outputData) return loadingDiv;

  return (
    <div className='output-wrapper'>
      <div>
        <p>Select which output data to view:</p>
        <Dropdown
          list={outputData.map((o) => o.name)}
          selected_index={0}
          onSelect={handleOnSelectOutput}
        ></Dropdown>
      </div>
      <div>
        <p>Select a view option:</p>
        <Dropdown
          list={Object.keys(viewOptions)}
          selected_index={0}
          onSelect={handleOnSelectView}
        ></Dropdown>
      </div>
      <div>
        {outputDatum ? 'Some data here' : 'No data yet for that output'}
      </div>
    </div>
  );
}
