import React, {useState, useCallback} from 'react';

import Grid from '@material-ui/core/Grid';
import Card from '@material-ui/core/Card';
import CardHeader from '@material-ui/core/CardHeader';
import CardContent from '@material-ui/core/CardContent';
import CardActions from '@material-ui/core/CardActions';
import Typography from '@material-ui/core/Typography';
import Button from '@material-ui/core/Button';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';

import {
  listDirectory,
  enqueueDirectoryFiles
} from '../../api/polyphemus/ingest_api';

import HostSelector from './ingest-dialog/host-selector';
import DirectorySelector from './ingest-dialog/directory-selector';
import FilesList, {IngestFile} from './ingest-dialog/files-list';

const IngestToBucketDialog = () => {
  const [files, setFiles] = useState([] as IngestFile[]);
  const [selectedFileNames, setSelectedFileNames] = useState([] as string[]);
  const [fetchError, setFetchError] = useState('');
  const [enqueueError, setEnqueueError] = useState('');
  const [host, setHost] = useState('');
  const [directory, setDirectory] = useState('');
  const invoke = useActionInvoker();

  const fetchFiles = useCallback(() => {
    setFetchError('');
    listDirectory(host, directory)
      .then((data: {[key: string]: IngestFile[]}) => {
        setFiles(data.files);
        setSelectedFileNames(data.files.map((f) => f.name));
      })
      .catch((e: Promise<string>) => {
        e.then((err: string) => setFetchError(err));
      });
  }, [host, directory]);

  const enqueueFiles = useCallback(() => {
    setEnqueueError('');
    enqueueDirectoryFiles(host, directory, selectedFileNames)
      .then(() => {
        invoke({type: 'DISMISS_DIALOG'});
      })
      .catch((e: Promise<string>) => {
        e.then((err: string) => setEnqueueError(err));
      });
  }, [host, directory, selectedFileNames]);

  const onChangeHost = useCallback(
    (newValue) => {
      setHost(newValue);
      if (files.length > 0) setFiles([]);
    },
    [files]
  );

  const onChangeDirectory = useCallback(
    (newValue) => {
      setDirectory(newValue);
      if (files.length > 0) setFiles([]);
    },
    [files]
  );

  return (
    <Grid
      container
      direction='column'
      className='ingest-dialog'
      alignItems='stretch'
    >
      <Grid item>
        <Card>
          <CardContent>
            <Grid container direction='column'>
              <Grid item>
                <HostSelector onChangeHost={onChangeHost} />
              </Grid>
              <Grid item>
                <DirectorySelector
                  host={host}
                  onChangeDirectory={onChangeDirectory}
                />
              </Grid>
              {fetchError && (
                <Grid item>
                  <Typography color='error'>{fetchError}</Typography>
                </Grid>
              )}
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
      <Grid item className='ingest-files-list-container'>
        <Card className='ingest-files-card'>
          <CardHeader
            subheader={`${files.length} file${
              1 === files.length ? '' : 's'
            } found`}
          />
          <CardContent className='ingest-files-list'>
            <FilesList
              allFiles={files}
              selectedFileNames={selectedFileNames}
              onChangeFileNames={setSelectedFileNames}
            />
            {enqueueError && (
              <Grid item>
                <Typography color='error'>{enqueueError}</Typography>
              </Grid>
            )}
          </CardContent>
          <CardActions>
            <Button
              variant='contained'
              color='primary'
              onClick={() => enqueueFiles()}
              disabled={
                !host || !directory || 0 === files.length || '' !== fetchError
              }
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
