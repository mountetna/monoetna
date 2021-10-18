import React, { useState, useEffect, useCallback, useContext } from 'react';
import { json_get } from 'etna-js/utils/fetch';
import { getDocuments } from 'etna-js/api/magma_api';

import Typography from '@material-ui/core/Typography';
import {makeStyles} from '@material-ui/core/styles';
import Grid from '@material-ui/core/Grid';
import Button from '@material-ui/core/Button';
import EtlConfig from './etl/etl-config';
import EtlCreate from './etl/etl-create';

import { MagmaContext } from 'etna-js/contexts/magma-context';

const useStyles = makeStyles((theme) => ({
  etls: {
    minWidth: '800px'
  },
  title: {
    padding: '10px'
  }
}));

const PolyphemusMain = ({project_name}:{project_name: string}) => {
  const [ etls, setEtls ] = useState< Etl[] >([]);
  const [ jobs, setJobs ] = useState< Job[] | null>(null);
  const [ create, setCreate ] = useState(false);
  const { models, setModels } = useContext(MagmaContext);

  const addEtl = (etl:Etl) => {
    let index = etls.findIndex( e => e.name == etl.name );
    let new_etls = (index == -1) ? etls.concat(etl) :
      etl.archived ? etls.filter( (e, i) => i != index ) :
      etls.map( (e,i) => i == index ? etl : e );

    setEtls(new_etls);
  }

  const classes = useStyles();

  useEffect( () => {
    json_get(`/api/etl/${project_name}/configs`).then(setEtls)
    json_get(`/api/etl/jobs`).then(setJobs);
    getDocuments({
      project_name,
      model_name: 'all',
      record_names: [],
      attribute_names: 'all',
    }, fetch).then(
      ({models}) => setModels(models)
    ).catch( e => console.log({e}) )
  }, [] );


  return <Grid id='polyphemus-main'>
    {
      !jobs ? null :
        <Grid className={classes.etls} item xs={12}>
          <Typography className={classes.title} variant="h5">
            {project_name} Data Loaders
          </Typography>
          {
            etls.map( (etl:Etl) => <EtlConfig key={etl.name} { ...etl } onUpdate={ addEtl } job={ jobs.find( j => j.name == etl.etl ) } />)
          }
          <Grid>
            <Button onClick={ () => setCreate(true) }>Add Loader</Button>
            <EtlCreate
              project_name={project_name}
              open={create}
              onClose={() => setCreate(false)}
              onCreate={addEtl}
              jobs={ jobs }/>
          </Grid>
        </Grid>
    }
  </Grid>;
}

export default PolyphemusMain;
