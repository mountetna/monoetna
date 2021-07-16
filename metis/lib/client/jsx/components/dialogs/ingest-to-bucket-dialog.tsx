import React, {useState, useCallback} from 'react';

import Grid from '@material-ui/core/Grid';
import Card from '@material-ui/core/Card';
import CardHeader from '@material-ui/core/CardHeader';
import CardContent from '@material-ui/core/CardContent';
import CardActions from '@material-ui/core/CardActions';
import Typography from '@material-ui/core/Typography';
import TextField from '@material-ui/core/TextField';
import Button from '@material-ui/core/Button';
import DescriptionIcon from '@material-ui/icons/Description';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';

import {
  listDirectory,
  queueDirectoryFiles
} from '../../api/polyphemus/ingest_api';

type IngestFile = {
  name: string;
  host: string;
};

const IngestToBucketDialog = ({}) => {
  const [files, setFiles] = useState([] as IngestFile[]);
  const [host, setHost] = useState('');
  const [directory, setDirectory] = useState('');
  const invoke = useActionInvoker();

  const fetchFiles = useCallback(() => {
    listDirectory(host, directory).then((data) => {
      setFiles(data.files);
    });
  }, [host, directory]);

  const queueFiles = useCallback(() => {
    queueDirectoryFiles(host, directory).then(() => {
      invoke({type: 'DISMISS_DIALOG'});
    });
  }, [host, directory]);

  return (
    <Grid container direction='column' className='ingest-dialog'>
      <Grid item>
        <Card>
          <CardContent>
            <Grid container direction='column'>
              <Grid item>
                <TextField
                  required
                  fullWidth
                  id='ingest-host'
                  label='Ingest host'
                  onChange={(e) => setHost(e.target.value)}
                />
              </Grid>
              <Grid item>
                <TextField
                  required
                  fullWidth
                  id='ingest-directory'
                  label='Directory path'
                  onChange={(e) => setDirectory(e.target.value)}
                />
              </Grid>
            </Grid>
          </CardContent>
          <CardActions>
            <Button
              variant='contained'
              color='primary'
              onClick={() => fetchFiles()}
              disabled={!host || !directory}
            >
              List files
            </Button>
          </CardActions>
        </Card>
      </Grid>
      <Grid item>
        <Card>
          <CardHeader
            subheader={`${files.length} file${
              1 === files.length ? '' : 's'
            } found`}
          />
          <CardContent className='ingest-files-list'>
            {files.map((file: IngestFile) => (
              <Typography>
                <DescriptionIcon />
                {file.name}
              </Typography>
            ))}
          </CardContent>
          <CardActions>
            <Button
              variant='contained'
              color='primary'
              onClick={() => queueFiles()}
              disabled={!host || !directory || 0 === files.length}
            >
              Queue files for ingestion!
            </Button>
          </CardActions>
        </Card>
      </Grid>
    </Grid>
  );
};

export default IngestToBucketDialog;
