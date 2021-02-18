import React, {useState, useContext, useEffect, useRef} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from 'etna-js/actions/message_actions';
import {addConsignment} from 'etna-js/plots/actions/manifest_actions';

import Dropdown from 'etna-js/components/inputs/dropdown';

import {VulcanContext} from '../../contexts/vulcan';
import {getData} from '../../api/vulcan';

import Plot from './output_renderers/plot';
import Raw from './output_renderers/raw';
import {stepDataUrls, stepOutputs} from '../../selectors/workflow_selector';

const loadingDiv = (
  <div className='browser'>
    <div id='loader-container'>
      <div className='loader'>Loading...</div>
    </div>
  </div>
);

// These are stubs that should be fetched from the context, normally.
// const WORKFLOW = require('../../../../server/data/steps.json');
// const UMAP_DATA = require('../../../../server/data/umap_data.json');

export default function Output() {
  const outputRef = useRef(null);
  const invoke = useActionInvoker();
  const {workflow, pathIndex, stepIndex, setData, setWorkflow} = useContext(
    VulcanContext
  );

  if (!workflow) return null;

  const [outputData, setOutputData] = useState(null);
  const [outputIndex, setOutputIndex] = useState(null);
  const [viewIndex, setViewIndex] = useState(null);
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

  // Used for DEV only
  // useEffect(() => {
  //   const STUB_URL = `https://vulcan.development.local/api/${CONFIG.project_name}/workflows/umap/umap_data`;

  //   WORKFLOW.steps[0][5].data_url = STUB_URL;

  //   setWorkflow(WORKFLOW);
  //   setOutputData(UMAP_DATA);
  // }, []);

  // useEffect(() => {
  //   // When the workflow updates, let's pull out
  //   //   the data for the current step.
  //   setOutputData(stepOutputs(workflow, 0, 4));
  // }, []);

  const viewOptions = {
    plot: Plot,
    raw: Raw
  };

  useEffect(() => {
    if (outputData) {
      let dataKey = Object.keys(outputData)[outputIndex];
      if (dataKey) {
        setOutputDatum(outputData[dataKey]);
      }
    }
  }, [outputData, outputIndex]);

  useEffect(() => {
    if (viewIndex) {
      // We should change this hash to be the actual cell hash.
      invoke(addConsignment('some-cell-hash', outputData));
    }
  }, [viewIndex]);

  let Component = viewOptions[Object.keys(viewOptions)[viewIndex]] || null;

  function handleOnSelectView(index) {
    setViewIndex(index);
  }

  function handleOnSelectOutput(index) {
    setOutputIndex(index);
  }

  if (!outputData) return loadingDiv;

  return (
    <div className='step-output'>
      <div>
        <p>Select a view option:</p>
        <Dropdown
          list={Object.keys(viewOptions)}
          selected_index={viewIndex}
          onSelect={handleOnSelectView}
        ></Dropdown>
      </div>
      <div ref={outputRef}>
        {Component ? (
          <Component md5sum='some-cell-hash' parentRef={outputRef}></Component>
        ) : null}
      </div>
    </div>
  );
}
