import React, { useState, useEffect, useCallback, useRef } from 'react';
import { json_get } from 'etna-js/utils/fetch';

import Typography from '@material-ui/core/Typography';
import {makeStyles} from '@material-ui/core/styles';
import Grid from '@material-ui/core/Grid';
import EtlConfig from './polyphemus-etl';

const useStyles = makeStyles((theme) => ({
  title: {
    padding: '10px'
  }
}));

const schemaFor = (etl: string, jobs: Job[]) : (any | null) => {
  let job = jobs.find( j => j.name == etl )
  return job ? job.schema : null;
}

const PolyphemusMain = ({project_name}:{project_name: string}) => {
  const [ etls, setEtls ] = useState< Etl[] | null>(null);
  const [ jobs, setJobs ] = useState< Job[] | null>(null);

  const classes = useStyles();

  useEffect( () => {
    json_get(`/api/${project_name}/etl/configs`).then(setEtls)
    json_get(`/api/${project_name}/etl/jobs`).then(setJobs)
  }, [] );
  return <Grid id='polyphemus-main'>
    {
      !(jobs && etls) ? null :
        <Grid item xs={6}>
          <Typography className={classes.title} variant="h5">
            {project_name} Data Loaders
          </Typography>
          {
            etls.map( (etl:Etl) => <EtlConfig key={etl.name} { ...etl } status='completed' schema={schemaFor(etl.etl, jobs)} />)
          }
        </Grid>
    }
  </Grid>;
}

export default PolyphemusMain;
