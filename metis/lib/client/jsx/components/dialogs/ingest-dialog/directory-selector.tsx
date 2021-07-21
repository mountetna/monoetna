import React, {useContext} from 'react';

import TextField from '@material-ui/core/TextField';
import Autocomplete from '@material-ui/lab/Autocomplete';

import {IngestContext} from '../../../contexts/ingest_context';

const DirectorySelector = ({
  host,
  onChangeDirectory
}: {
  host: string;
  onChangeDirectory: (value: string) => void;
}) => {
  const {state} = useContext(IngestContext);

  if (!host || !state.hosts[host]) return null;

  return (
    <Autocomplete
      id='ingest-directory'
      options={state.hosts[host].directories.sort()}
      getOptionLabel={(option: string) => option}
      fullWidth
      renderInput={(params: any) => (
        <TextField {...params} label='Ingest directory' variant='outlined' />
      )}
      classes={{
        popper: 'ingest-directory-options',
        root: 'ingest-selector'
      }}
      onChange={(e: React.ChangeEvent<{}>, value: string | null) =>
        onChangeDirectory(value || '')
      }
    />
  );
};

export default DirectorySelector;
