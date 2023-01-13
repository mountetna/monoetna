import React, { useState, createContext, useContext, useCallback, useEffect } from 'react';
import {ModalDialogContainer} from 'etna-js/components/ModalDialogContainer';
import {metisPath} from 'etna-js/api/metis_api.js'
import Autocomplete from '@material-ui/lab/Autocomplete';
import { TextField } from '@material-ui/core';
import { json_get } from 'etna-js/utils/fetch';

function dispValue(value: string | null) {
  return value == null ? '' : value;
}

export function pickBucket({ project_name="example", setBucket, initialValue }: {
  project_name?: string;
  setBucket: any;
  initialValue?: string;
}){
  const [bucketList, setBucketList] = useState([] as string[]);
  const [bucket, setBucketInternal] = useState(initialValue? initialValue : null);
  const [inputState, setInputState] = useState(dispValue(bucket));

  useEffect( () => {
    json_get(metisPath(`${project_name}/list`)).then(
      bucket_list => setBucketList(bucket_list)
    );
  }, [] );

  function onChangeAction(event: any, e: string | null) {
    setBucketInternal(e)
    setBucket(e)
  }
  
  console.log({bucketList})
  console.log({bucket})

  return (
    <Autocomplete
      clearOnBlur={true}
      options={bucketList}
      value={bucket}
      onChange={onChangeAction}
      inputValue={inputState}
      onInputChange={(event: any, newInputState: string) => {
        setInputState(newInputState);
      }}
      style={{minWidth: 100, paddingTop: 8}}
      renderInput={(params) => (
        <TextField
          {...params}
          error={ bucket==null || inputState != dispValue(bucket)}
          label="Pick a Bucket"
          size='small'
          InputLabelProps={{shrink: true}}
        />
      )}
    />
  );
}


