import React, { useState, useCallback, useEffect, useMemo } from 'react';
import {metisPath} from 'etna-js/api/metis_api.js'
import Autocomplete from '@material-ui/lab/Autocomplete';
import { Grid, TextField } from '@material-ui/core';
import { json_get } from 'etna-js/utils/fetch';
import SubdirectoryArrowRightIcon from '@material-ui/icons/SubdirectoryArrowRight';
import FolderSharpIcon from '@material-ui/icons/FolderSharp';
import {useAsyncCallback} from 'etna-js/utils/cancellable_helpers';

declare const CONFIG: {[key: string]: any};

declare type folderSummary = {folders: string[], files: string[]}

function nullToEmptyString(value: string | null | undefined) {
  return value == null ? '' : value;
}

function addSlash(path: string) {
  if (path.slice(-1) != '/') {
    path = path + '/'
  }
  return path
}
function stripSlash(path: string) {
  if (path.slice(-1) == '/') {
    path = path.slice(0,-1)
  }
  return path
}

function arrayToPath(pathElements: string[], basePath: string) {
  const base = basePath!='' ? addSlash(basePath) : ''
  return base + pathElements.join('/')
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

async function getFolderContents(path: string, project_name: string, bucket: string) {
  const full_path = `${bucket}/${path}`
  let contents = {} as folderSummary
  await json_get(metisPath(`${project_name}/list/${full_path}`)).then(
    metis_return => {
      contents = simplifyFolderListReturn(metis_return)
    }
  )
  return contents
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
    <Grid container direction='row'>
      <Grid item container alignItems='flex-end' style={{width: 'auto'}}>
        {'\uD83E\uDEA3'}
      </Grid>
      <Grid item style={{flex: '1 1 auto'}}>
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
      </Grid>
    </Grid>
  );
}

function SingleFileOrFolder({ optionSet, path, setTarget, onEmpty, basePath, topLevelPlaceholer, topLevelLabel, autoOpen, awaitingContents}: {
  optionSet: string[];
  path: string;
  setTarget: Function;
  onEmpty: Function;
  autoOpen: boolean;
  basePath: string;
  topLevelPlaceholer?: string;
  topLevelLabel: string;
  awaitingContents: boolean;
}){
  const [target, setTargetInternal] = useState('');
  const [inputState, setInputState] = useState(target);

  return <Autocomplete
    key={path+'-selection'}
    disablePortal
    openOnFocus={autoOpen}
    autoHighlight
    options={optionSet.concat('')}
    loading={awaitingContents}
    inputValue={inputState}
    value={target}
    onInputChange={(event: any, newInputState: string) => {
      setInputState(newInputState);
    }}
    onChange={ (event: any, e: string | null) => {
      const nextTarget = nullToEmptyString(e)
      setTargetInternal(nextTarget)
      if (nextTarget == '') {
        onEmpty()
      } else {
        setTarget(nextTarget)
      }
    }}
    filterOptions={(options, state) => {
      const query = inputState != target ? inputState : '';
      return optionSet.filter((o) => {
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
}) {
  const [pathArray, setPathArray] = useState(pathToArray(path, basePath));
  const labelUse = label!==undefined ? label : allowFiles ? "File or Folder" : "Folder"
  const [contentsSeen, setContentsSeen] = useState({} as {[k: string]: folderSummary})
  const [awaitingContents, setAwaitingContents] = useState(true)

  const contentUse = (folderContents: folderSummary) => {
    return allowFiles ?
      folderContents.folders.concat(folderContents.files) :
      folderContents.folders
  }
  const contentKey = useCallback( (folderPath) => {
    return `${project_name}_${bucket}_${folderPath}`
  }, [project_name, bucket])

  const pathSeen = useCallback( (folderPath) => {
    return Object.keys(contentsSeen).includes(contentKey(folderPath))
  }, [contentKey, contentsSeen])
  const subPathsSeen = useCallback( (folderPath) => {
    const pathKey = contentKey(folderPath)
    if (!pathSeen(folderPath)) return false
    if (contentsSeen[pathKey].folders.length == 0) {
      return true
    } else {
      return contentsSeen[pathKey].folders.every(
        (folder: string) => {
          pathSeen(`${folderPath}/${folder}`)
        }
      )
    }
  }, [contentKey, contentsSeen])
  
  const addNewConstentsSeen = (folderPath: string, newContents: folderSummary) => {
    const fullContentsSeen = {...contentsSeen}
    const key = contentKey(folderPath)
    fullContentsSeen[key] = newContents
    setContentsSeen(fullContentsSeen)
  }
  const [fetchFolderContents] = useAsyncCallback(function* (
    path: string, callback: Function
  ) {
    const full_path = `${bucket}/${path}`
    const metis_return = yield json_get(metisPath(`${project_name}/list/${full_path}`));
    callback(metis_return);
  },
  [project_name, bucket, contentsSeen]);

  // Reset and give time to re-buffer contents if bucket is changed
  useEffect(() => {
    setAwaitingContents(true)
    setPathArray(pathToArray('', basePath))
  }, [bucket])

  // Buffer current folder contents
  useEffect(() => {
    const newPath = stripSlash(arrayToPath(pathArray, basePath))
    const newKey = contentKey(newPath)
    // Obtain contents of current folder
    if (!pathSeen(newPath)) {
      fetchFolderContents(newPath, (x: any) => {
        const contents = simplifyFolderListReturn(x)
        addNewConstentsSeen(newPath, contents)
        setAwaitingContents(false)
      })
    } else {
      setAwaitingContents(false)
    }
  }, [pathArray])

  // Buffer potentially next folder contents
  useEffect( () => {
    const newPath = stripSlash(arrayToPath(pathArray, basePath))
    if (pathSeen(newPath) && !subPathsSeen(newPath)) {
      const newKey = contentKey(newPath)
      contentsSeen[newKey].folders.forEach(
        (folder: string) => {
          const subPath = newPath=='' ? folder : `${newPath}/${folder}`
          if (!pathSeen(subPath)) {
            fetchFolderContents(subPath, (x: any) => {
              addNewConstentsSeen(subPath, simplifyFolderListReturn(x))
            })
          }
        }
      )
    }
  }, [pathArray, contentsSeen])

  // Update main readout (path) whenever anything updates pathArray
  useEffect(() => {
    setPath(arrayToPath(pathArray, basePath));
  }, [pathArray])

  // onChange for inners
  const updateTarget = (newPath: string, currentPathSet: string[], depth: number) => {
    // Also trim at the level just picked in case not actually the deepest
    let nextPathSet = [...currentPathSet].slice(0,depth+1)
    nextPathSet[depth] = newPath
    if (contentUse(contentsSeen[contentKey(arrayToPath(nextPathSet, basePath))]).length > 0 ) {
      nextPathSet.push('')
    }
    setAwaitingContents(true)
    setPathArray(nextPathSet)
  }
  
  // clear for inners
  const trimPath = (depth: number) => {
    const nextPathSet = depth > 0 ? [...pathArray].slice(0,depth) : ['']
    setPathArray(nextPathSet)
  }

  // console.log({pathArray})
  // console.log({contentsSeen})

  return !pathSeen(stripSlash(arrayToPath(pathArray, basePath))) ? null : <Grid 
    key={`file-or-folder-picker-${project_name}-${bucket}`}
    container
    direction='column'
  >
    {pathArray.map( (p, index, fullSet) => {
      const before = index > 1 ? pathArray[index-1] : '';
      const icon = index > 0 ? <SubdirectoryArrowRightIcon fontSize='small'/> : '\uD83D\uDCC1';
      const thisPath = arrayToPath(fullSet.slice(0,index), basePath)
      const ready = (index+1) == fullSet.length ? !awaitingContents : true
      return (
        <Grid item container
          key={`sub-picker-${project_name}-${bucket}-${before}`+ index}
          >
          <Grid item container alignItems='flex-end' style={{width: 'auto'}}>
            {icon}
          </Grid>
          <Grid item style={{flex: '1 1 auto'}}>
            <SingleFileOrFolder
              optionSet={ ready ? contentUse(contentsSeen[contentKey(thisPath)]) : []}
              awaitingContents={ ready }
              setTarget={(e: string) => updateTarget(e, pathArray, index)}
              onEmpty={ () => { trimPath(index) } }
              path={thisPath}
              autoOpen={index!=0}
              basePath={basePath}
              topLevelPlaceholer={topLevelPlaceholer}
              topLevelLabel={labelUse}
            />
          </Grid>
        </Grid>
      )
    })}
  </Grid>
}


