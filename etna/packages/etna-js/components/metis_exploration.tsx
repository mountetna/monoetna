import React, { useState, useCallback, useEffect, useMemo } from 'react';
import {metisPath} from 'etna-js/api/metis_api.js'
import Autocomplete from '@material-ui/lab/Autocomplete';
import { Grid, TextField } from '@material-ui/core';
import { json_get } from 'etna-js/utils/fetch';
import SubdirectoryArrowRightIcon from '@material-ui/icons/SubdirectoryArrowRight';
import FolderSharpIcon from '@material-ui/icons/FolderSharp';

declare const CONFIG: {[key: string]: any};

function nullToEmptyString(value: string | null | undefined) {
  return value == null ? '' : value;
}

function addSlash(path: string) {
  if (path.slice(-1) != '/') {
    path = path + '/'
  }
  return path
}

function arrayToPath(pathElements: string[], basePath: string) {
  return addSlash(basePath) + pathElements.join('/')
}

function pathToArray(path: string, basePath: string) {
  return path.replace(addSlash(basePath), '').split("/")
}

function samePath(path1: string, path2: string) {
  return addSlash(path1) == addSlash(path2)
}

function simplifyFolderListReturn(fullList: any) {
  return {
    folders: fullList.folders.map( (f: any) => f.folder_name) as string[],
    files: fullList.files.map( (f: any) => f.file_name) as string[]
  }
}

export function PickBucket({ project_name=CONFIG.project_name, setBucket, bucket, label="Bucket"}: {
  project_name?: string;
  setBucket: any;
  bucket: string;
  label?: string;
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
      style={{paddingTop: 8}}
      renderInput={(params) => (
        <TextField
          {...params}
          error={ bucket==null || bucket == '' || inputState != nullToEmptyString(bucket)}
          label={label}
          size='small'
          InputLabelProps={{shrink: true}}
        />
      )}
    />
  );
}

function SingleFileOrFolder({ project_name, bucket, path, setTarget, onEmpty, basePath, topLevelPlaceholer, topLevelLabel, allowFiles = true, autoOpen}: {
  project_name?: string;
  bucket: string;
  path: string;
  setTarget: Function;
  onEmpty: Function;
  allowFiles?: boolean;
  autoOpen: boolean;
  basePath: string;
  topLevelPlaceholer?: string;
  topLevelLabel: string;
}){
  const [fullTargetList, setFullTargetList] = useState({folders: [] as string[], files: [] as string[]})
  function restrictTargetList(fullList = fullTargetList, allow_files = allowFiles) {
    return allow_files ? fullList.folders.concat(fullList.files) : fullList.folders
  }
  const [targetsWithContents, setTargetsWithContents] = useState([] as string[])
  const [targetList, setTargetList] = useState(restrictTargetList(fullTargetList));
  const [target, setTargetInternal] = useState('');
  const [inputState, setInputState] = useState(target);

  useEffect( () => {
    if (bucket != null) {
      const full_path = path.length>0 ? `${bucket}/${path}` : bucket
      json_get(metisPath(`${project_name}/list/${full_path}`)).then(
        metis_return => {
          const fullList = simplifyFolderListReturn(metis_return)
          setFullTargetList(fullList)
          setTargetList(restrictTargetList(fullList))
          let withContents = [] as string[]
          fullList.folders.forEach(
            (folder) => {
              json_get(metisPath(`${project_name}/list/${full_path}/${folder}`)).then(
                child_metis_return => {
                  const contents_allowed = allowFiles ?
                    child_metis_return.folders.concat(child_metis_return.files) :
                    child_metis_return.folders
                  if (contents_allowed.length > 0 ) {
                    withContents.push(folder)
                  }
              })
          })
          setTargetsWithContents(withContents)
        }
      );
    }
    // console.log(path)
    // console.log(targetList)
    // console.log(targetsWithContents)
  }, [bucket, project_name, path] );

  return project_name == null ? null : (
    <Autocomplete
      key={path+'-selection'}
      disablePortal
      openOnFocus={autoOpen}
      autoHighlight
      options={targetList.concat('')}
      inputValue={inputState}
      value={target}
      onInputChange={(event: any, newInputState: string) => {
        setInputState(newInputState);
      }}
      onChange={ (event: any, e: string | null) => {
        const nextTarget = nullToEmptyString(e)
        setTargetInternal(nextTarget)
        setTarget(nextTarget.concat(targetsWithContents.includes(nextTarget) ? '/' : ''))
        if (nextTarget == '') {
          onEmpty()
        }
      }}
      filterOptions={(options, state) => {
        const query = inputState != target ? inputState : '';
        return targetList.filter((o) => {
          return query == null ? true : o.indexOf(query) > -1;
        });
      }}
      style={{paddingTop: samePath(path, basePath) ? 8 : undefined}}
      renderInput={(params) => (
        <TextField
          {...params}
          error={ inputState != nullToEmptyString(target) }
          autoFocus
          placeholder={samePath(path, basePath) ? topLevelPlaceholer : undefined}
          label={samePath(path, basePath) ? topLevelLabel : undefined}
          size='small'
          InputLabelProps={{shrink: true}}
        />
      )}
    />
  );
}

export function PickFileOrFolder({ project_name=CONFIG.project_name, bucket, setPath, path, allowFiles, label, basePath, topLevelPlaceholer=''}: {
  project_name?: string;
  bucket: string;
  setPath: Function;
  path: string; // the overall "value" / output
  label?: string;
  basePath: string; // an immutable portion of path.  Can be '' to access the entire bucket and to be compatible with exploring across buckets in sync with a PickBucket companion.
  topLevelPlaceholer?: string;
  allowFiles?: boolean
}){
  const [pathArray, setPathArray] = useState(pathToArray(path, basePath));
  const labelUse = label!==undefined ? label : allowFiles ? "File or Folder" : "Folder"

  useEffect(() => {
    setPathArray(pathToArray(path, basePath));
  }, [path])

  // onChange for inners
  const updateTarget = useCallback(
    (newPath: string, currentPathSet: string[], depth: number) => {
      // Also trim at the level just picked in case not actually the deepest
      // console.log(newPath)
      let nextPathSet = [...currentPathSet].slice(0,depth+1)
      nextPathSet[depth] = newPath
      setPath(arrayToPath(nextPathSet, basePath))
    }, [setPath])
  
  // clear for inners
  const trimPath = useCallback(
    (depth: number) => {
    const nextPathSet = depth > 0 ? [...pathArray].slice(0,depth+1) : ['']
    setPath(arrayToPath(nextPathSet, basePath))
  }, [setPath])

  const targetSelectors = useMemo(() => {
    return pathArray.map( (p, index, fullSet) => {
      const before = index > 1 ? pathArray[index-1] : '';
      const icon = index > 0 ? <SubdirectoryArrowRightIcon fontSize='small'/> : <FolderSharpIcon fontSize='small'/>;
      return (
        <Grid item container
          key={`sub-picker-${bucket}-${project_name}-${before}`+ index}
          >
          <Grid item container alignItems='flex-end' style={{width: 'auto'}}>
            {icon}
          </Grid>
          <Grid item style={{flex: '1 1 auto'}}>
            <SingleFileOrFolder
              project_name={project_name}
              bucket={bucket}
              setTarget={(e: string) => updateTarget(e, pathArray, index)}
              onEmpty={ () => { trimPath(index) } }
              path={arrayToPath(fullSet.slice(0,index), basePath)}
              allowFiles={allowFiles}
              autoOpen={index!=0}
              basePath={basePath}
              topLevelPlaceholer={topLevelPlaceholer}
              topLevelLabel={labelUse}
            />
          </Grid>
        </Grid>
      )
    })
  }, [pathArray, project_name, bucket, updateTarget])

  // console.log(path)
  // console.log(pathArray)

  return bucket == null ? null : <Grid 
    key={`file-or-folder-picker-${bucket}-${project_name}`}
    container
    direction='column'
  >
    {targetSelectors}
  </Grid>
}


