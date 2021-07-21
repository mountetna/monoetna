import React, {useContext} from 'react';

import TextField from '@material-ui/core/TextField';
import Autocomplete from '@material-ui/lab/Autocomplete';

import {IngestContext, Host} from '../../../contexts/ingest_context';

const HostSelector = ({
  onChangeHost
}: {
  onChangeHost: (host: string) => void;
}) => {
  const {state} = useContext(IngestContext);

  if (!state.hosts) return null;

  function prettify(host: Host) {
    return `(${host.alias}) ${host.host}`;
  }

  return (
    <Autocomplete
      id='ingest-host'
      options={Object.values(state.hosts).sort()}
      getOptionLabel={(option: Host) => prettify(option)}
      fullWidth
      renderInput={(params: any) => (
        <TextField {...params} label='Ingest host' variant='outlined' />
      )}
      classes={{
        popper: 'ingest-host-options',
        root: 'ingest-selector'
      }}
      onChange={(e: React.ChangeEvent<{}>, value: Host | null) =>
        onChangeHost(value ? value.host : '')
      }
    />
  );
};

export default HostSelector;
