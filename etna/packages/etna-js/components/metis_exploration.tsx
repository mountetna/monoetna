import React, { useState, createContext, useContext, useCallback, useEffect, useMemo } from 'react';
import {metisPath} from 'etna-js/api/metis_api.js'
import Autocomplete from '@material-ui/lab/Autocomplete';
import { Button, makeStyles, TextField } from '@material-ui/core';
import { json_get } from 'etna-js/utils/fetch';
import Tooltip from '@material-ui/core/Tooltip';
import DeleteIcon from '@material-ui/icons/Delete';
import AddIcon from '@material-ui/icons/Add';
import SubdirectoryArrowRightIcon from '@material-ui/icons/SubdirectoryArrowRight';

const useStyles = makeStyles((theme) => ({
  button: {
    margin: '0.5rem'
  }
}));

declare const CONFIG: {[key: string]: any};

function nullToEmptyString(value: string | null | undefined) {
  return value == null ? '' : value;
}

function arrayToPath(pathElements: string[]) {
  return pathElements.join('/')
}

function pathToArray(path: string) {
  return path.split("/")
}

export function PickBucket({ project_name=CONFIG.project_name, setBucket, bucket}: {
  project_name?: string;
  setBucket: any;
  bucket: string;
}){
  const [bucketList, setBucketList] = useState([] as string[]);
  const [inputState, setInputState] = useState(nullToEmptyString(bucket));

  useEffect( () => {
    json_get(metisPath(`${project_name}/list`)).then(
      metis_return => {setBucketList(metis_return.buckets.map( (b: any) => b.bucket_name))}
    );
  }, [project_name] );

  return (
    <Autocomplete
      options={bucketList}
      value={nullToEmptyString(bucket)}
      onChange={ (event: any, e: string | null) => {
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

function SingleFileOrFolder({ project_name, bucket, path, setTarget, allowFiles = true, allowFolders = true }: {
  project_name?: string;
  bucket: string;
  path: string;
  setTarget: any;
  allowFiles?: boolean; 
  allowFolders?: boolean;
}){
  const [fullTargetList, setFullTargetList] = useState({folders: [] as any, files: [] as any})
  function restrictTargetList(fullList = fullTargetList, allow_files = allowFiles, allow_folders = allowFolders) {
    let out = allow_folders ? fullList.folders.map(f => f.folder_name) : []
    return allow_files ? out.concat(fullList.files.map(f => f.file_name)) : out
  }
  const [targetList, setTargetList] = useState([''].concat(restrictTargetList(fullTargetList)));
  const [target, setTargetInternal] = useState('');
  const [inputState, setInputState] = useState(target);

  useEffect( () => {
    if (bucket != null) {
      const full_path = path.length>0 ? `${bucket}/${path}` : bucket
      json_get(metisPath(`${project_name}/list/${full_path}`)).then(
        metis_return => {
          setFullTargetList(metis_return)
          setTargetList([''].concat(restrictTargetList(metis_return)))
        }
      );
    }
    // console.log({path})
    // console.log({targetList})
  }, [bucket, project_name, path] );

  return project_name == null ? null : (
    <Autocomplete
      key={path+'-selection'}
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

export function PickFileOrFolder({ project_name=CONFIG.project_name, bucket, setPath, path }: {
  project_name?: string;
  bucket: string;
  setPath: Function;
  path: string;
}){
  // const originalTarget = nullToEmptyString(initialValue);
  // const originalPath = pathToArray(originalTarget);
  const [pathArray, setPathArray] = useState(pathToArray(nullToEmptyString(path)));

  useEffect(() => {
    setPathArray(pathToArray(nullToEmptyString(path)));
  }, [path])
  const classes = useStyles();

  // onChange for inners
  const updateTarget = useCallback(
    (newPath: string, currentPathSet: string[], depth: number) => {
      // Also trim at the level just picked in case not actually the deepest
      let nextPathSet = [...currentPathSet].slice(0,depth+1)
      nextPathSet[depth] = newPath
      setPath(arrayToPath(nextPathSet))
    }, [setPath])
  
  // clear for inners
  const trimPath = useCallback(
    (depth: number) => {
    const nextPathSet = depth > 0 ? [...pathArray].slice(0,depth) : ['']
    setPath(arrayToPath(nextPathSet))
  }, [setPath])

  const targetSelectors = useMemo(() => {
    return pathArray.map( (p, index, fullSet) => {
      const indent = index > 1 ? (index-1)*10 : 0;
      const icon = index > 0 ? <SubdirectoryArrowRightIcon/> : null;
      return (
        <div 
          style={{display: 'inline-flex', paddingLeft: indent}}
          key={`sub-picker-${bucket}-${project_name}-`+ index}
          >
          {icon}
          <SingleFileOrFolder
            project_name={project_name}
            bucket={bucket}
            setTarget={(e: string) => updateTarget(e, pathArray, index)}
            path={arrayToPath(fullSet.slice(0,index))}
          />
          <Tooltip title='Clear or Remove Level'>
            <Button
              startIcon={<DeleteIcon />}
              onClick={ () => { trimPath(index)} }
              size='small'
              color='primary'
              className={classes.button}
            >
            </Button>
          </Tooltip>
        </div>
      )
    })
  }, [pathArray, project_name, bucket, updateTarget])

  // console.log({pathSet: pathArray})
  // console.log({bucket})
  // console.log({project_name})

  return bucket == null ? null : <div key={`file-or-folder-picker-${bucket}-${project_name}`}>
    {targetSelectors}
    {/* Show button only if end is a folder with internal contents & at least one folder is fileAllow=false */}
    <Tooltip title='Add Level'>
      <Button
        startIcon={<AddIcon />}
        onClick={() => { setPath(path+'/') }}
        size='small'
        color='primary'
        className={classes.button}
      >
      </Button>
    </Tooltip>
  </div>
}


