import React, { useState, createContext, useContext, useCallback, useEffect } from 'react';
import {metisPath} from 'etna-js/api/metis_api.js'
import Autocomplete from '@material-ui/lab/Autocomplete';
import { TextField } from '@material-ui/core';
import { json_get } from 'etna-js/utils/fetch';

declare const CONFIG: {[key: string]: any};

function dispValue(value: string | null) {
  return value == null ? '' : value;
}

export function PickBucket({ project_name, setBucket, initialValue }: {
  project_name?: string;
  setBucket: any;
  initialValue?: string;
}){
  const [bucketList, setBucketList] = useState([] as string[]);
  const [bucket, setBucketInternal] = useState(initialValue? initialValue : null);
  const [inputState, setInputState] = useState(dispValue(bucket));

  useEffect( () => {
    if (project_name == null) {
      project_name = CONFIG.project_name
    };
    json_get(metisPath(`${project_name}/list`)).then(
      metis_return => {setBucketList(metis_return.buckets.map( (b: any) => b.bucket_name))}
    );
  }, [project_name] );
  
  console.log({bucketList})
  console.log({bucket})

  return (
    <Autocomplete
      options={bucketList}
      value={dispValue(bucket)}
      onChange={ (event: any, e: string | null) => {
        setBucketInternal(e)
        setBucket(e)
      }}
      disableClearable
      disablePortal
      inputValue={inputState}
      onInputChange={(event: any, newInputState: string) => {
        setInputState(newInputState);
      }}
      style={{minWidth: 100, paddingTop: 8}}
      renderInput={(params) => (
        <TextField
          {...params}
          error={ bucket==null || bucket == '' || inputState != dispValue(bucket)}
          label="Bucket"
          size='small'
          InputLabelProps={{shrink: true}}
        />
      )}
    />
  );
}

