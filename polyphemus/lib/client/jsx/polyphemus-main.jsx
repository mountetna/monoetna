import React, { useState, useEffect, useCallback } from 'react';
import { json_get } from 'etna-js/utils/fetch';

const PolyphemusMain = ({project_name}) => {
  const [ etls, setEtls ] = useState(null);

  useEffect( () => {
    json_get(`/api/${project_name}/etls`).then((es) => {
      setEtls(es)
    })
  }, [] );
  return <div id='polyphemus-main'>
    {project_name}
    { etls ? <div className='etls'>
      {
        etls.map( etl => <div className='etl' key={etl.etl}>
          { etl.name }
        </div> )
      }
    </div> : null }
  </div>;
}

export default PolyphemusMain;
