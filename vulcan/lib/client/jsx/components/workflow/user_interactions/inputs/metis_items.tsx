import React, {useEffect} from 'react';
import {WithInputParams} from './input_types';
import {some, withDefault} from '../../../../selectors/maybe';
import { PickFileOrFolder, PickBucket } from 'etna-js/components/metis_exploration';
import {useSetsDefault} from './useSetsDefault';
import { Grid } from '@material-ui/core';

declare const CONFIG: {[key: string]: any};
type metisPathType = { bucket: string, path: string, type: 'file' | 'folder' | null }

function _MetisLocationInput({onChange, label, allowFiles, data, ...props}: WithInputParams<{
    label?: string,
    allowFiles: boolean
    }, metisPathType, metisPathType
  >) {
  const value = useSetsDefault({bucket: '', path: '', type: null}, props.value, onChange);

  function updateKeyToVal(key: 'bucket' | 'path' | 'type', val: string | null, fullValues = {...value}) {
    // console.log('setting ', key, ' to ', val)
    const newValue = fullValues
    newValue[key] = val as any
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

  console.log({value})

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
          disablePortal={false}
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
          allowFiles={allowFiles}
          label={label?  label+', '+'File/Folder' : undefined}
          basePath={''}
          topLevelPlaceholder='Path to target file'
          disablePortal={false}
        />
      </Grid>
    </Grid>
    
  );
}

export function MetisFileInput({onChange, label, data, ...props}: WithInputParams<{label?: string}, metisPathType, metisPathType>) {
  return _MetisLocationInput({onChange, label, allowFiles: true, data, ...props})
}
export function MetisFolderInput({onChange, label, data, ...props}: WithInputParams<{label?: string}, metisPathType, metisPathType>) {
  return _MetisLocationInput({onChange, label, allowFiles: false, data, ...props})
}