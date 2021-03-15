import React from 'react';

import Icon from 'etna-js/components/icon';

// To get webpack to pick up the files.
require('../../../img/umap.png');
require('../../../img/default.png');
require('../../../img/test_workflow.png');

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
        <img src={`/images/${workflow.image || 'default.png'}`} />
      </figure>
      <div className='workflow-card-description'>
        <div className='row'>
          <div className='label'>Name</div>
          <div className='value'>{workflow.name}</div>
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
            <Icon className='card-tags' icon='tags'></Icon> Tags
          </div>
          <div className='value'>{MultiLineOutput(workflow.tags)}</div>
        </div>
      </div>
    </div>
  );
}
