import React, {useState, useCallback} from 'react';

import Grid from '@material-ui/core/Grid';
import Card from '@material-ui/core/Card';
import CardHeader from '@material-ui/core/CardHeader';
import CardContent from '@material-ui/core/CardContent';
import CardActions from '@material-ui/core/CardActions';
import Typography from '@material-ui/core/Typography';
import TextField from '@material-ui/core/TextField';
import Button from '@material-ui/core/Button';

import {listDirectory} from '../../api/polyphemus/ingest_api';

type IngestFile = {
  name: string;
  host: string;
};

const IngestToBucketDialog = ({}) => {
  const [files, setFiles] = useState([] as IngestFile[]);
  const [host, setHost] = useState('');
  const [directory, setDirectory] = useState('');

  const fetchFiles = useCallback(() => {
    listDirectory(host, directory).then((data) => {
      setFiles(data.files);
    });
  }, [host, directory]);

  return (
    <Grid container direction='column'>
      <Grid item>
        <Card>
          <CardContent>
            <TextField
              required
              id='ingest-host'
              label='Ingest host'
              onChange={(e) => setHost(e.target.value)}
            />
            <TextField
              required
              id='ingest-directory'
              label='Directory path'
              onChange={(e) => setDirectory(e.target.value)}
            />
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
          <CardHeader>
            <Typography>{files.length} files found</Typography>
          </CardHeader>
          <CardContent>
            {files.map((file: IngestFile) => {
              <Typography>{file.name}</Typography>;
            })}
          </CardContent>
        </Card>
      </Grid>
    </Grid>
  );
};

export default IngestToBucketDialog;
