import React, { useState, useEffect, useCallback, useRef } from 'react';
import { json_get } from 'etna-js/utils/fetch';

import Typography from '@material-ui/core/Typography';
import Card from '@material-ui/core/Card';
import CardContent from '@material-ui/core/CardContent';
import {makeStyles} from '@material-ui/core/styles';
import CardActions from '@material-ui/core/CardActions';
import Button from '@material-ui/core/Button';
import Box from '@material-ui/core/Box';
import Collapse from '@material-ui/core/Collapse';

import ConfigScript from './config-script';

const Config = ({config, schema}:{config: any, schema: any}) =>  {
  return <ConfigScript script={config} schema={schema}/>
}

type Job = {
  name: string;
  schema: any
};

type Etl = {
  name: string,
  etl: string,
  config: any
};

const schemaFor = (etl: string, jobs: Job[]) : (any | null) => {
  let job = jobs.find( j => j.name == etl )
  return job ? job.schema : null;
}

const PolyphemusMain = ({project_name}:{project_name: string}) => {
  const [ etls, setEtls ] = useState< Etl[] | null>(null);
  const [ jobs, setJobs ] = useState< Job[] | null>(null);
  const [ expanded, setExpanded ] = useState(false);

  const toggleExpanded = () => setExpanded(!expanded);

  useEffect( () => {
    json_get(`/api/${project_name}/etl/configs`).then(setEtls)
    json_get(`/api/${project_name}/etl/jobs`).then(setJobs)
  }, [] );
  return <div id='polyphemus-main'>
    {project_name}
    {
      !(jobs && etls) ? null :
        <div className='etls'>
          <span className='title'>ETLs</span>
          {
            etls.map( ({etl,name,config} : Etl) => <Card className='etl' key={etl}>
              <CardContent>
                <Typography component="h5" variant="h5">
                  { name }
                </Typography>

                <CardActions>
                  <Button variant="text" onClick={toggleExpanded} size="small">Configure</Button>
                </CardActions>
                <Collapse in={expanded} timeout="auto" unmountOnExit>
                  <CardContent>
                    <Config config={ config } schema={ schemaFor(etl, jobs) } />
                  </CardContent>
                </Collapse>
              </CardContent>
            </Card> )
          }
        </div>
    }
  </div>;
}

export default PolyphemusMain;
