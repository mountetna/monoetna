import React, { useState, useCallback, useEffect, useMemo } from 'react';
import {metisPath} from 'etna-js/api/metis_api.js'
import Autocomplete from '@material-ui/lab/Autocomplete';
import { Grid, IconButton, TextField } from '@material-ui/core';
import { json_get } from 'etna-js/utils/fetch';
import SubdirectoryArrowRightIcon from '@material-ui/icons/SubdirectoryArrowRight';
import AddIcon from '@material-ui/icons/Add';

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
  const base_use = addSlash(basePath)
  if (path.startsWith(base_use)) path = path.slice(base_use.length)
  return path.split("/")
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
    <Grid container direction='row'>
      <Grid item container alignItems='flex-end' style={{width: 'auto'}}>
        <i className="fa fa-trash" aria-hidden="true" style={{padding:'3px 5px 10px 3px'}}/>
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

function FileOrFolderInner({ optionSet, path, target, setTarget, onEmpty, placeholder, label, autoOpen}: {
  optionSet: string[];
  target: string;
  path: string;
  setTarget: Function;
  onEmpty: Function;
  autoOpen: boolean;
  placeholder?: string;
  label?: string;
}){
  const [inputState, setInputState] = useState(target);

  return <Autocomplete
    key={path+'-selection'}
    disablePortal
    openOnFocus={true}
    autoHighlight
    options={optionSet.concat('')}
    inputValue={inputState}
    value={target}
    onInputChange={(event: any, newInputState: string) => {
      setInputState(newInputState);
    }}
    onChange={ (event: any, e: string | null) => {
      const nextTarget = nullToEmptyString(e)
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
    style={{paddingTop: label!==undefined ? 8 : undefined}}
    renderInput={(params) => (
      <TextField
        {...params}
        error={ inputState != nullToEmptyString(target) }
        autoFocus={autoOpen}
        placeholder={placeholder}
        label={label}
        size='small'
        InputLabelProps={{shrink: true}}
      />
    )}
  />
}

export function PickFileOrFolder({ project_name=CONFIG.project_name, bucket, setPath, path, allowFiles=true, label, basePath, topLevelPlaceholder=''}: {
  project_name?: string;
  bucket: string;
  setPath: Function;
  path: string; // the overall "value" / output
  label?: string;
  basePath: string; // an immutable portion of path.  Can be '' to access the entire bucket and to be compatible with exploring across buckets in sync with a PickBucket companion.
  topLevelPlaceholder?: string;
  allowFiles?: boolean
}) {
  const [pathArray, setPathArray] = useState(pathToArray(path, basePath));
  const labelUse = label!==undefined ? label : allowFiles ? "File or Folder" : "Folder"
  const [showAddButton, setShowAddButton] = useState(-1)
  const [contentsSeen, setContentsSeen] = useState({} as {[k: string]: folderSummary | undefined | string})
  const [firstRender, setFirstRender] = useState(true) // tracked so we can leave newly rendered dropdowns closed when not the user's focus
  const [fetchContents, setFetchContents] = useState([] as string[])
  const [fetching, setFetching] = useState(false)

  const contentUse = (folderContents: folderSummary) => {
    return allowFiles ?
      folderContents.folders.concat(folderContents.files) :
      folderContents.folders
  }
  const pathInMetis = useCallback( (folderPath: string) => {
    return `${project_name}/list/${bucket}/${folderPath}`
  }, [project_name, bucket])

  const pathSeen = useCallback( (folderPath: string) => {
    const key = pathInMetis(folderPath)
    return key in contentsSeen && !["undefined","string"].includes(typeof(contentsSeen[key]))  // todo: use string for error cases when trying to fetch folder contents
  }, [pathInMetis, contentsSeen])
  
  // fetch contents from metis. Caveat: only works one-at-a-time because contentsSeen gets stale.
  const fetchFolderContents = useCallback( (path: string, callback?: (folderSummary: folderSummary) => void) => {
    const key = pathInMetis(path)
    setContentsSeen({...contentsSeen, [key]: undefined})
    json_get(metisPath(key)).then((metis_return) => {
      const folderSummary = simplifyFolderListReturn(metis_return)
      const newContentsSeen = {...contentsSeen, [key]: folderSummary}
      setContentsSeen(newContentsSeen)
      if (callback) callback(folderSummary)
    })
  },
  [project_name, bucket, contentsSeen]);

  // Reset and give time to re-buffer contents if bucket is changed
  useEffect(() => {
    if (Object.keys(contentsSeen).length>0) { // Skip at initialization
      setFirstRender(true)
      setPathArray(pathToArray('', basePath))
    }
  }, [bucket])

  useEffect(() => {
    // Update main readout (path) whenever anything updates pathArray
    setPath(stripSlash(arrayToPath(pathArray, basePath)));
  }, [pathArray])

  // Grab folder contents at very start
  useEffect( () => {
    let contentsNeeded = ['']
    for (let i=1; i <= pathArray.length; i++) {
      const thisPath = stripSlash(arrayToPath(pathArray.slice(0,i), basePath))
      contentsNeeded.push(thisPath)
    }
    setFetchContents(contentsNeeded)
  }, [])

  // Get contents, one at a time.
  useEffect( ()=>{
    if (fetchContents.length>0) {
      if (pathSeen(fetchContents[0])) {
        setFetching(false)
        setFetchContents(fetchContents.slice(1))
      } else if (!fetching) {
        setFetching(true)
        if (fetchContents[0]==path) {
          fetchFolderContents(fetchContents[0], (f) => { if (contentUse(f).length > 0) setPathArray([...pathArray, '']) })
        } else {
          fetchFolderContents(fetchContents[0])
        }
      }
    }
  }, [fetchContents, contentsSeen, pathArray])

  // onChange for inners
  const updateTarget = (newPath: string, currentPathSet: string[], depth: number) => {
    setFirstRender(false)
    // Also trim at the level just picked in case not actually the deepest
    let nextPathSet = [...currentPathSet].slice(0,depth+1)
    nextPathSet[depth] = newPath
    const targetPath = arrayToPath(nextPathSet, basePath)
    if (!pathSeen(targetPath)) {
      fetchFolderContents(targetPath, (f) => {
        if (contentUse(f).length > 0) {
          setPathArray([...nextPathSet, ''])
        } else {
          setPathArray(nextPathSet)
        }
      }) // todo: eliminate the double-fetch caused here & by fetchContents route once pathArray is updated
    } else {
      if (contentUse(contentsSeen[pathInMetis(targetPath)] as folderSummary).length > 0) {
        setPathArray([...nextPathSet, ''])
      } else {
        setPathArray(nextPathSet)
      }
    }
    setShowAddButton(-1)
  }
  
  // clear for inners
  const trimPath = (depth: number) => {
    const nextPathSet = depth > 0 ? [...pathArray].slice(0,depth) : ['']
    setPathArray(nextPathSet)
    if (depth > 0) setShowAddButton(depth-1)
  }

  return <Grid 
    key={`file-or-folder-picker-${project_name}-${bucket}`}
    container
    direction='column'
  >
    {pathArray.map( (p, index, fullArray) => {
      // const before = index > 1 ? fullSet[index-1] : '';
      const icon = index > 0 ? <SubdirectoryArrowRightIcon fontSize='small'/> : <i className="fas fa-folder" aria-hidden="true" style={{padding:'2px'}}/>;
      const thisPath = arrayToPath(fullArray.slice(0,index), basePath)
      const ready = pathSeen(thisPath)
      return (
        <Grid item container
          key={`sub-picker-${project_name}-${bucket}-`+ index}
          >
          <Grid item container alignItems='flex-end' style={{width: 'auto', paddingBottom: '7px', paddingRight: '2px'}}>
            {icon}
          </Grid>
          <Grid item style={{flex: '1 1 auto'}}>
            {ready && <FileOrFolderInner
              optionSet={ contentUse(contentsSeen[pathInMetis(thisPath)] as folderSummary)}
              target={p}
              setTarget={(e: string) => updateTarget(e, pathArray, index)}
              onEmpty={ () => { trimPath(index) } }
              path={thisPath}
              autoOpen={!firstRender}
              placeholder={(index==0) ? topLevelPlaceholder : undefined}
              label={(index==0) ? labelUse : undefined}
            />}
          </Grid>
          {showAddButton==index ? <Grid item container alignItems='flex-end' style={{width: 'auto'}}>
            <IconButton
              aria-label="Re-add next level"
              size="small"
              onClick={() => {
                setPathArray([...fullArray, ''])
                setShowAddButton(-1)
                }}>
              <AddIcon/>
            </IconButton>
          </Grid> : null
          }
        </Grid>
      )
    })}
  </Grid>
}

export function PickFolder({ project_name=CONFIG.project_name, bucket, setPath, path, label, basePath, topLevelPlaceholder=''}: {
  project_name?: string;
  bucket: string;
  setPath: Function;
  path: string; // the overall "value" / output
  label?: string;
  basePath: string; // an immutable portion of path.  Can be '' to access the entire bucket and to be compatible with exploring across buckets in sync with a PickBucket companion.
  topLevelPlaceholder?: string;
}) {
  return <PickFileOrFolder
    project_name={project_name}
    bucket={bucket}
    setPath={setPath}
    path={path}
    label={label}
    basePath={basePath}
    topLevelPlaceholder={topLevelPlaceholder}
    allowFiles={false}
    />
}

