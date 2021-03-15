import React from 'react';

import Icon from 'etna-js/components/icon';

// To get webpack to pick up the files.
require('../../../img/umap.png');
require('../../../img/default.png');
require('../../../img/test_workflow.png');

export default function Card({workflow}) {
  const displayFields = ['name', 'authors', 'lastModified'];
  return (
    <div className='workflow-card'>
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
          <div className='value'>
            {workflow.authors.map((a, ind) => (
              <div key={ind}>{a}</div>
            ))}
          </div>
        </div>
        <div className='row'>
          <div className='label'>Last Modified</div>
          <div className='value'>{workflow.lastModified}</div>
        </div>
        <div className='row'>
          <div className='label'>
            <Icon className='card-tags' icon='tags'></Icon> Tags
          </div>
          <div className='value'>
            {workflow.tags.map((t, ind) => (
              <div key={ind}>{t}</div>
            ))}
          </div>
        </div>
      </div>
    </div>
  );
}
