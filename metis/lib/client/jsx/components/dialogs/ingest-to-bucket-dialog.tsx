import React, { useState, useCallback, useContext } from 'react';

import Grid from '@material-ui/core/Grid';
import Card from '@material-ui/core/Card';
import CardHeader from '@material-ui/core/CardHeader';
import CardContent from '@material-ui/core/CardContent';
import CardActions from '@material-ui/core/CardActions';
import Typography from '@material-ui/core/Typography';
import TextField from '@material-ui/core/TextField';
import Button from '@material-ui/core/Button';
import DescriptionIcon from '@material-ui/icons/Description';
import Autocomplete from '@material-ui/lab/Autocomplete';

import { useActionInvoker } from 'etna-js/hooks/useActionInvoker';

import { IngestContext, Host } from '../../contexts/ingest_context';

import {
  listDirectory,
  enqueueDirectoryFiles
} from '../../api/polyphemus/ingest_api';

type IngestFile = {
  name: string;
  host: string;
};

const HostSelector = ({ onChangeHost }: { onChangeHost: (host: string) => void }) => {

  const { state } = useContext(IngestContext);

  if (!state.hosts) return null;

  function prettify(host: Host) {
    return `(${host.alias}) ${host.host}`;
  }

  return <Autocomplete
    id='ingest-host'
    options={Object.values(state.hosts).sort()}
    getOptionLabel={(option: Host) => prettify(option)}
    fullWidth
    renderInput={(params: any) => (
      <TextField
        {...params}
        label='Ingest host'
        variant='outlined'
      />
    )}
    classes={{
      popper: 'ingest-host-options',
      root: 'ingest-selector'
    }}
    onChange={(e: React.ChangeEvent<{}>, value: Host | null) =>
      onChangeHost(value ? value.host : '')
    }
  />
}

const DirectorySelector = ({ host, onChangeDirectory }: { host: string, onChangeDirectory: (value: string) => void }) => {
  const { state } = useContext(IngestContext);

  if (!host || !state.hosts[host]) return null;

  return <Autocomplete
    id='ingest-directory'
    options={state.hosts[host].directories.sort()}
    getOptionLabel={(option: string) => option}
    fullWidth
    renderInput={(params: any) => (
      <TextField
        {...params}
        label='Ingest directory'
        variant='outlined'
      />
    )}
    classes={{
      popper: 'ingest-directory-options',
      root: 'ingest-selector'
    }}
    onChange={(e: React.ChangeEvent<{}>, value: string | null) =>
      onChangeDirectory(value || '')
    }
  />
}


const IngestToBucketDialog = () => {
  const [files, setFiles] = useState([] as IngestFile[]);
  const [fetchError, setFetchError] = useState('');
  const [enqueueError, setEnqueueError] = useState('');
  const [host, setHost] = useState('');
  const [directory, setDirectory] = useState('');
  const invoke = useActionInvoker();
  const { state } = useContext(IngestContext);

  const fetchFiles = useCallback(() => {
    setFetchError('');
    listDirectory(host, directory)
      .then((data: { [key: string]: IngestFile[] }) => {
        setFiles(data.files);
      })
      .catch((e: Promise<string>) => {
        e.then((err: string) => setFetchError(err));
      });
  }, [host, directory]);

  const enqueueFiles = useCallback(() => {
    setEnqueueError('');
    enqueueDirectoryFiles(host, directory)
      .then(() => {
        invoke({ type: 'DISMISS_DIALOG' });
      })
      .catch((e: Promise<string>) => {
        e.then((err: string) => setEnqueueError(err));
      });
  }, [host, directory]);

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
    <Grid container direction='column' className='ingest-dialog' alignItems='stretch'>
      <Grid item>
        <Card>
          <CardContent>
            <Grid container direction='column'>
              <Grid item>
                <HostSelector onChangeHost={onChangeHost} />
              </Grid>
              <Grid item>
                <DirectorySelector host={host} onChangeDirectory={onChangeDirectory} />
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
            subheader={`${files.length} file${1 === files.length ? '' : 's'
              } found`}
          />
          <CardContent className='ingest-files-list'>
            {files.map((file: IngestFile) => (
              <Typography>
                <DescriptionIcon />
                {file.name}
              </Typography>
            ))}
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
