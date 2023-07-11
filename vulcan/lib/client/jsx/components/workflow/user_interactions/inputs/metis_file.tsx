import React, {useEffect} from 'react';
import {WithInputParams} from './input_types';
import {some, withDefault} from '../../../../selectors/maybe';
import { PickFileOrFolder, PickBucket } from 'etna-js/components/metis_exploration';
import {useSetsDefault} from './useSetsDefault';
import {selectDefaultBoolean} from './monoids';
import { Grid } from '@material-ui/core';
import { Checkbox, FormControlLabel, FormControlLabelProps } from '@material-ui/core';
import { metisPath } from 'etna-js/api/metis_api.js'
import { json_get } from 'etna-js/utils/fetch';

declare const CONFIG: {[key: string]: any};
type metisPathType = { bucket: string, path: string, type: 'file' | 'folder' | null }

export default function FileInput({onChange, label, data, ...props}: WithInputParams<{label?: string}, metisPathType, metisPathType>) {
  // console.log({props})
  const value = useSetsDefault({bucket: '', path: '', type: null}, props.value, onChange);
  // console.log({value})

  function updateKeyToVal(key: string, val: string | null, fullValues = {...value}) {
    // console.log('setting ', key, ' to ', val)
    const newValue = fullValues
    newValue[key] = val
    onChange(some(newValue))
  }

  function updateType(t: 'file' | 'folder' | null) {
    updateKeyToVal('type', t)
  }

  useEffect( () => {
    if (value.bucket == '' && value.type != null) {
      updateType(null)
    } 
  }, [value.bucket])

  return (
    !value || value.bucket == undefined ? null :
    <Grid container direction='column'>
      <Grid item>
        <PickBucket
          setBucket={ (b: string) => {
            updateKeyToVal('bucket', b)
          }}
          bucket={value.bucket}
          label={label? label+', '+'Bucket' : undefined}
        />
      </Grid>
      <Grid>
        <PickFileOrFolder
          bucket={value.bucket}
          setPath={ (p: string) => {
            updateKeyToVal('path', p)
          }}
          path={value.path}
          useTargetType={updateType}
          target
          allowFiles={true}
          label={label?  label+', '+'File/Folder' : undefined}
          basePath={''}
          topLevelPlaceholder='Path to target file'
        />
      </Grid>
    </Grid>
    
  );
}