import React from 'react';

import Icon from 'etna-js/components/icon';

import {workflowName} from '../../selectors/workflow_selectors';
import ImageMemo from './image_memo';

// To get webpack to pick up the files.
require('../../../img/umap.png');
require('../../../img/barplot.png');
require('../../../img/scatter.png');
require('../../../img/yplot.png');
require('../../../img/default.png');
require('../../../img/add_integers.png');

function MultiLineOutput(iterator) {
  return iterator.map((item, ind) => (
    <div key={ind} className='sub_item'>
      {item}
    </div>
  ));
}

export default function Card({workflow, onClick}) {
  return (
    <div className='workflow-card' onClick={onClick}>
      <figure className='workflow-card-image'>
        <ImageMemo
          src={`/images/${workflow.image || 'default.png'}`}
          alt='Workflow image'
        />
      </figure>
      <div className='workflow-card-description'>
        <div className='row'>
          <div className='label'>Name</div>
          <div className='value'>
            {workflow.displayName || workflowName(workflow)}
          </div>
        </div>
        <div className='row'>
          <div className='label'>Authors</div>
          <div className='value'>{MultiLineOutput(workflow.authors)}</div>
        </div>
        <div className='row'>
          <div className='label'>Last Modified</div>
          <div className='value'>{workflow.lastModified}</div>
        </div>
        <div className='row'>
          <div className='label'>
            <Icon className='card-tags' icon='tags'/> Tags
          </div>
          <div className='value'>{MultiLineOutput(workflow.tags)}</div>
        </div>
      </div>
    </div>
  );
}
