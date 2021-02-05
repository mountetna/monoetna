import React, {useState, useContext, useEffect} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from 'etna-js/actions/message_actions';

import {VulcanContext} from '../../contexts/vulcan';
import {getData} from '../../api/vulcan';

import Plot from './output_renderers/plot';
import Raw from './output_renderers/raw';
import {stepDataUrls} from '../../selectors/workflow_selector';

export default function Output() {
  const invoke = useActionInvoker();
  const {workflow, pathIndex, stepIndex, setData} = useContext(VulcanContext);

  useEffect(() => {
    let dataUrls = stepDataUrls(workflow, pathIndex);
    dataUrls.forEach((url) => {
      getData(url)
        .then((data) => {
          setData(url, data);
        })
        .catch((e) => {
          invoke(showMessages(e));
        });
    });
  }, [pathIndex, stepIndex]);

  const viewOptions = {
    plot: Plot,
    raw: Raw
  };

  let Component;
  function handleOnSelect(index) {
    Component = viewOptions[Object.keys(viewOptions)[index]];
  }

  return (
    <div className='output-wrapper'>
      <Dropdown
        list={Object.keys(viewOptions)}
        selected_index={0}
        onSelect={handleOnSelect}
      ></Dropdown>
      <div></div>
    </div>
  );
}
