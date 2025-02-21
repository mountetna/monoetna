import React, {
  useState,
  useEffect,
  useCallback,
  useContext,
  useMemo
} from 'react';
import Grid from '@material-ui/core/Grid';
import TextField from '@material-ui/core/TextField';

import {makeStyles, Theme} from '@material-ui/core/styles';
import {SchemaContext, SchemaProvider} from './schema-context';

import {PickBucket} from 'etna-js/components/metis_exploration';
import {Script, Job} from '../polyphemus';
import AddModel from './add-model';

export const useStyles = makeStyles((theme: Theme) => ({
}));

export type Config = {
  bucket_name: string,
  magic_string: string,
  deposit_root_path: string,
  metis_root_path: string,
  ingest_root_path: string
};

const IngestionForm = ({
  config,
  project_name,
  job,
  update
}: {
  project_name: string;
  config: Config;
  job: Job | undefined;
  update: (config: Config) => void;
}) => {
  const classes = useStyles();

  const {setSchema} = useContext(SchemaContext);

  useEffect(() => {
    if (job) setSchema(job.schema);
  }, [job]);

  const bucket_name = config.bucket_name || '';

  const { magic_string, metis_root_path, ingest_root_path, deposit_root_path } = config;

  const setBucket = useCallback(
    (bucket_name: string) => {
      let c = {...config, bucket_name };
      update(c as Config);
    },
    [config, update]
  )

  return (
    <Grid container className={classes.form} direction='column'>
      <Grid item container>
        <Grid item xs={3} style={{paddingLeft:'10px'}}>Magic string</Grid>
        <Grid item xs={9}>
          <TextField
            fullWidth
            placeholder='Match string for target files'
            onChange={(e) => update({ ...config, magic_string: e.target.value})}
            value={magic_string}
          />
        </Grid>
      </Grid>
      <Grid item container>
        <Grid item xs={3} style={{paddingLeft:'10px'}}>Ingest root path</Grid>
        <Grid item xs={9}>
          <TextField
            fullWidth
            placeholder='Search path for files on ingest host'
            onChange={(e) => update({ ...config, ingest_root_path: e.target.value})}
            value={ingest_root_path}
          />
        </Grid>
      </Grid>
      <Grid item container>
        <Grid item xs={3} style={{paddingLeft:'10px'}}>Deposit root path</Grid>
        <Grid item xs={9}>
          <TextField
            fullWidth
            placeholder='Path to copy files onto deposit host'
            onChange={(e) => update({ ...config, deposit_root_path: e.target.value})}
            value={deposit_root_path}
          />
        </Grid>
      </Grid>
      <Grid item container>
        <Grid item xs={3} style={{paddingLeft:'10px'}}>Metis bucket</Grid>
        <Grid item xs={9}>
          <PickBucket setBucket={setBucket} project_name={project_name} bucket={bucket_name} label={null}/>
        </Grid>
      </Grid>
      <Grid item container>
        <Grid item xs={3} style={{paddingLeft:'10px'}}>Metis root path</Grid>
        <Grid item xs={9}>
          <TextField
            fullWidth
            placeholder='Path to copy files onto Metis'
            onChange={(e) => update({ ...config, metis_root_path: e.target.value})}
            value={metis_root_path}
          />
        </Grid>
      </Grid>
    </Grid>
  );
};

const IngestionFormProvider = (props: any) => (
  <SchemaProvider>
    <IngestionForm {...props} />
  </SchemaProvider>
);

export default IngestionFormProvider;
