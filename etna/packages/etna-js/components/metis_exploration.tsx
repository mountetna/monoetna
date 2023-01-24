import React, { useState, createContext, useContext, useCallback, useEffect, useMemo } from 'react';
import {metisPath} from 'etna-js/api/metis_api.js'
import Autocomplete from '@material-ui/lab/Autocomplete';
import { Button, TextField } from '@material-ui/core';
import { json_get } from 'etna-js/utils/fetch';

declare const CONFIG: {[key: string]: any};

function nullToEmptyString(value: string | null) {
  return value == null ? '' : value;
}

function arrayToPath(pathElements: string[]) {
  return pathElements.join('/')
}

export function PickBucket({ project_name=CONFIG.project_name, setBucket, initialValue }: {
  project_name?: string;
  setBucket: any;
  initialValue?: string;
}){
  const [bucketList, setBucketList] = useState([] as string[]);
  const [bucket, setBucketInternal] = useState(initialValue? initialValue : null);
  const [inputState, setInputState] = useState(nullToEmptyString(bucket));

  useEffect( () => {
    json_get(metisPath(`${project_name}/list`)).then(
      metis_return => {setBucketList(metis_return.buckets.map( (b: any) => b.bucket_name))}
    );
  }, [project_name] );
  
  // console.log({bucketList})
  // console.log({bucket})

  return (
    <Autocomplete
      options={bucketList}
      value={nullToEmptyString(bucket)}
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
          error={ bucket==null || bucket == '' || inputState != nullToEmptyString(bucket)}
          label="Bucket"
          size='small'
          InputLabelProps={{shrink: true}}
        />
      )}
    />
  );
}

function SingleFileOrFolder({ project_name, bucket, path, setTarget, key, initialValue, allowFiles = true, allowFolders = true }: {
  project_name?: string;
  bucket: string;
  path: string;
  setTarget: any;
  key: number;
  initialValue?: string;
  allowFiles?: boolean; 
  allowFolders?: boolean;
}){
  const [fullTargetList, setFullTargetList] = useState({folders: [] as any, files: [] as any})
  function restrictTargetList(fullList = fullTargetList, allow_files = allowFiles, allow_folders = allowFolders) {
    let out = allow_folders ? fullList.folders.map(f => f.folder_name) : []
    return allow_files ? out.concat(fullList.files.map(f => f.file_name)) : out
  }
  const [targetList, setTargetList] = useState(restrictTargetList(fullTargetList));
  const [target, setTargetInternal] = useState(initialValue? initialValue : '');
  const [inputState, setInputState] = useState(nullToEmptyString(target));

  useEffect( () => {
    if (bucket != null) {
      const full_path = path.length>0 ? `${bucket}/${path}` : bucket
      json_get(metisPath(`${project_name}/list/${full_path}`)).then(
        metis_return => {
          setFullTargetList(metis_return)
          setTargetList(restrictTargetList(metis_return))
        }
      );
    }
  }, [bucket, project_name, path] );
  
  console.log({fullTargetList})
  console.log({targetList})
  console.log({target})

  return project_name == null ? null : (
    <Autocomplete
      key={`sub-picker-${bucket}-${project_name}-`+ key}
      options={targetList}
      value={target}
      onChange={ (event: any, e: string | null) => {
        setTargetInternal(nullToEmptyString(e))
        setTarget(nullToEmptyString(e))
      }}
      disableClearable
      disablePortal
      inputValue={inputState}
      onInputChange={(event: any, newInputState: string) => {
        setInputState(newInputState);
      }}
      // style={{minWidth: 100}}
      renderInput={(params) => (
        <TextField
          {...params}
          error={ target==null || inputState != nullToEmptyString(target)}
          size='small'
          InputLabelProps={{shrink: true}}
        />
      )}
    />
  );
}

export function PickFileOrFolder({ project_name=CONFIG.project_name, bucket, setTarget, initialValue }: {
  project_name?: string;
  bucket: string;
  setTarget: any;
  initialValue?: string;
}){
  const defaultPath = initialValue? [initialValue] : ['']
  const [pathSet, setPathSet] = useState(defaultPath);

  // Clear if the project or bucket change
  useEffect( () => {
    setPathSet(defaultPath)
  }, [bucket, project_name] );

  // onChange for inner levels
  const updateTarget = useCallback(
    (newPath: string, currentPathSet: string[], depth: number) => {
      // Also trim at the level just picked in case not actually the deepest
      let nextPathSet = [...currentPathSet].slice(0,depth+1)
      nextPathSet[depth] = newPath
      setPathSet([...nextPathSet])
      setTarget(arrayToPath(nextPathSet))
    }, [setTarget])

  const targetSelectors = useMemo(() => {
    return pathSet.map( (p, index) => <SingleFileOrFolder
      project_name={project_name}
      bucket={bucket}
      initialValue={p}
      key={index}
      setTarget={(e: string) => updateTarget(e, pathSet, index)}
      path={arrayToPath(pathSet.slice(0,index))}
    />)
  }, [pathSet, project_name, bucket, updateTarget])
  let displaySet = [] as any[]
  targetSelectors.forEach((e,i) => {
    displaySet = displaySet.concat(i==0 ? [e] : ['/', e])
  })

  console.log({pathSet})
  console.log({bucket})
  console.log({project_name})

  return bucket == null ? null : <div>
    {displaySet.map( x => x)}
    {/* Show button only if end is a folder with internal contents & at least one folder is fileAllow=false */}
    <Button
      onClick={() => { setPathSet(pathSet.concat([''])) }}
      >
      +
    </Button>
  </div>
}


