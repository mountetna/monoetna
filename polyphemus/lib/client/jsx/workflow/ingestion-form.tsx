import React, {
  useState,
  useEffect,
  useCallback,
  useContext,
  useMemo
} from 'react';
import Grid from '@material-ui/core/Grid';
import Tabs from '@material-ui/core/Tabs';
import Tab from '@material-ui/core/Tab';
import Typography from '@material-ui/core/Typography';
import Tooltip from '@material-ui/core/Tooltip';
import Select from '@material-ui/core/Select';
import MenuItem from '@material-ui/core/MenuItem';
import IconButton from '@material-ui/core/IconButton';
import Button from '@material-ui/core/Button';
import AddIcon from '@material-ui/icons/Add';
import DeleteIcon from '@material-ui/icons/Delete';
import CopyIcon from '@material-ui/icons/FileCopy';
import Pagination from '@material-ui/lab/Pagination';
import Switch from '@material-ui/core/Switch';
import FormControlLabel from '@material-ui/core/FormControlLabel';

import {makeStyles, Theme} from '@material-ui/core/styles';
import {SchemaContext, SchemaProvider} from './schema-context';

import {Script, Job} from '../polyphemus';
import AddModel from './add-model';

export const useStyles = makeStyles((theme: Theme) => ({
}));

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
        regex
      </Grid>
      <Grid item container>
        root_dir
      </Grid>
      <Grid item container>
        file_regex
      </Grid>
      <Grid item container>
        sftp_root_dir
      </Grid>
      <Grid item container>
        path_to_write_files
      </Grid>
      <Grid item container>
        bucket_name
      </Grid>
      <Grid item container>
        metis_root_path
      </Grid>
      <Grid item container>
        c4_root_path
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
