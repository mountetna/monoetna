import React, { useState, useCallback, useEffect, useMemo } from 'react';
import {metisPath} from 'etna-js/api/metis_api.js'
import Autocomplete from '@material-ui/lab/Autocomplete';
import Grid from '@material-ui/core/Grid';
import IconButton from '@material-ui/core/IconButton';
import TextField from '@material-ui/core/TextField';
import Tooltip from '@material-ui/core/Tooltip';
import AutorenewIcon from '@material-ui/icons/Autorenew';
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

function simplifyFolderListReturn(fullList: any) {
  return {
    folders: fullList.folders.map( (f: any) => f.folder_name) as string[],
    files: fullList.files.map( (f: any) => f.file_name) as string[]
  }
}

const contentUse = (folderContents: folderSummary, allowFiles: boolean) => {
  return allowFiles ?
    folderContents.folders.concat(folderContents.files) :
    folderContents.folders
}

export function PickBucket({ project_name=CONFIG.project_name, setBucket, bucket, label="Bucket", className, disablePortal=true}: {
  project_name?: string;
  setBucket: any;
  bucket: string;
  label?: string | null;
  className?: string;
  disablePortal?: boolean;
}){
  const [bucketList, setBucketList] = useState([] as string[]);
  const [inputState, setInputState] = useState(nullToEmptyString(bucket));

  useEffect( () => {
    json_get(metisPath(`${project_name}/list`)).then(
      metis_return => {setBucketList(metis_return.buckets.map( (b: any) => b.bucket_name))}
    );
  }, [project_name] );

  return (
    <Grid container direction='row' className={className}>
      <Grid item container alignItems='flex-end' style={{width: 'auto'}}>
        <i className="fa fa-trash" aria-hidden="true" style={{padding:'5px 5px 5px 3px'}}/>
      </Grid>
      <Grid item style={{flex: '1 1 auto'}}>
        <Autocomplete
          // Add '' to avoid error message
          options={bucketList.concat('')}
          // But don't show '' as option
          filterOptions={(options, state) => {
            const query = inputState != bucket ? inputState : '';
            return bucketList.filter((o) => {
              return query == null ? true : o.indexOf(query) > -1;
            });
          }}
          value={nullToEmptyString(bucket)}
          onChange={ (event: any, e: string | null) => {
            setBucket(e)
          }}
          disableClearable
          disablePortal={disablePortal}
          inputValue={inputState}
          onInputChange={(event: any, newInputState: string) => {
            setInputState(newInputState);
          }}
          style={{paddingTop: label!=undefined ? 8 : undefined}}
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

function FileOrFolderInner({ optionSet, path, target, setTarget, onEmpty, placeholder, label, autoOpen, disablePortal=true}: {
  optionSet: string[];
  target: string;
  path: string;
  setTarget: Function;
  onEmpty: Function;
  autoOpen: boolean;
  placeholder?: string;
  label?: string;
  disablePortal?: boolean;
}){
  const [inputState, setInputState] = useState(target);

  return <Autocomplete
    key={path+'-selection'}
    disablePortal={disablePortal}
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

export function PickFileOrFolder({ project_name=CONFIG.project_name, bucket, setPath, path, useTargetType, allowFiles=true, label, basePath, topLevelPlaceholder='', className, disablePortal=true}: {
  project_name?: string;
  bucket: string;
  setPath: Function;
  path: string; // the overall "value" / output
  useTargetType?: (type: 'folder' | 'file' | null)=>void; // null = not yet determined OR not getting tracked
  label?: string;
  basePath: string; // an immutable portion of path.  Can be '' to access the entire bucket and to be compatible with exploring across buckets in sync with a PickBucket companion.
  topLevelPlaceholder?: string;
  allowFiles?: boolean;
  className?: string;
  disablePortal?: boolean;
}) {
  const [firstBucketPath, setFirstBucketPath] = useState([bucket,path]);
  const [pathArray, setPathArray] = useState(pathToArray(path, basePath));
  const [targetType, setTargetType] = useState(null as 'folder' | 'file' | null) // null = not yet determined
  const labelUse = label!==undefined ? label : allowFiles ? "File or Folder" : "Folder"
  const [showAddButton, setShowAddButton] = useState(-1) // holds index of pathArray level needing plus button shown next to it
  const [contentsSeen, setContentsSeen] = useState({} as {[k: string]: folderSummary | undefined})
  const [firstRender, setFirstRender] = useState(true) // tracked so we can leave newly rendered dropdowns closed when not the user's focus
  const [fetchContents, setFetchContents] = useState([] as string[])
  const [fetching, setFetching] = useState(false)
  const [reset, setReset] = useState(false);

  const pathInMetis = useCallback( (folderPath: string) => {
    return `${project_name}/list/${bucket}/${folderPath}`
  }, [project_name, bucket])
  const pathSeen = useCallback( (folderPath: string) => {
    const key = pathInMetis(folderPath)
    return key in contentsSeen && !["undefined","string"].includes(typeof(contentsSeen[key]))  // todo: use string for error cases when trying to fetch folder contents
  }, [pathInMetis, contentsSeen])
  
  function useType(type: 'file' | 'folder' | null, doTransmit: boolean = true) {
    if (doTransmit && useTargetType) useTargetType(type)
    if (type != targetType) setTargetType(type)
    return type
  }
  const useDetermineType = useCallback( (target: string, containerPath: string, doUse: boolean = true, doTransmit: boolean = true)=> {
    const containerMetisPath = pathInMetis(containerPath);
    const targetContainerContent = contentsSeen[containerMetisPath];
    let newType = null as 'file' | 'folder' | null;

    // console.log({contentsSeen})
    // console.log("Looking for ", target, " in ", containerMetisPath)
    // console.log({targetContainerContent})

    if (targetContainerContent == undefined) {
      if (!fetchContents.includes(containerPath)) {
        setFetchContents([...fetchContents, containerPath])
      }
      // newType = null;
    } else {
      const whichContents: ('folders'|'files')[] = (['folders', 'files'] as ('folders'|'files')[])
        .filter( contentType => targetContainerContent[contentType].includes(target) )
      newType = whichContents.length>0 ? whichContents[0].slice(0, -1) as 'file' | 'folder' : null
    }
    return doUse ? useType(newType, doTransmit) : newType
  }, [contentsSeen, fetchContents, pathInMetis])

  // Ensure type is determined even if containing folder contents weren't known at time of path update
  useEffect( () => {
    if (targetType == null && fetchContents.length == 0) {
      const containerPathArray = pathArray.slice(0,-1)
      const target = pathArray[pathArray.length-1]
      // without passing upwards to not trigger extraneous updates in parent component updates at page refresh
      if (target && target != '') {
        useDetermineType(target, stripSlash(arrayToPath(containerPathArray, basePath)), true, false)
      }
      if (target && target == '') {
        useType('folder', false)
      }
    }
  }, [targetType, fetchContents])
  
  // fetch contents from metis. Caveat: only works one-at-a-time because contentsSeen gets stale.
  const fetchFolderContents = useCallback( (path: string, callback?: (folderSummary: folderSummary) => void) => {
    const key = pathInMetis(path)
    // console.log("fetching contents of ", key)
    setContentsSeen({...contentsSeen, [key]: undefined})
    json_get(metisPath(key)).then((metis_return) => {
      const folderSummary = simplifyFolderListReturn(metis_return)
      const newContentsSeen = {...contentsSeen, [key]: folderSummary}
      setContentsSeen(newContentsSeen)
      if (callback) callback(folderSummary)
    })
  },
  [project_name, bucket, contentsSeen]);

  // Determine folder contents to grab at very start, or once first valid bucket given 
  useEffect( () => {
    if (bucket != '' && Object.keys(contentsSeen).length==0 && !fetching) {
      let contentsNeeded = ['']
      for (let i=1; i <= pathArray.length-1; i++) {
        const thisPath = stripSlash(arrayToPath(pathArray.slice(0,i), basePath))
        contentsNeeded.push(thisPath)
      }
      setFetchContents(contentsNeeded)
    }
  }, [bucket])

  // Get contents of paths in fetchContents, one at a time.
  useEffect( ()=>{
    if (fetchContents.length>0) {
      if (pathSeen(fetchContents[0])) {
        setFetching(false)
        setFetchContents(fetchContents.slice(1))
      } else if (!fetching) {
        setFetching(true)
        if (path !='' && fetchContents[0]==path) {
          fetchFolderContents(fetchContents[0], (f) => { if (contentUse(f, allowFiles).length > 0) setPathArray([...pathArray, '']) })
        } else {
          fetchFolderContents(fetchContents[0])
        }
      }
    }
  }, [fetchContents, contentsSeen, pathArray])

  // Reset if bucket is changed and don't auto-open dropdowns
  useEffect(() => {
    if (Object.keys(contentsSeen).length>0) { // Skip at initialization
      setFirstRender(true)
      setPathArray(pathToArray('', basePath))
      setFetchContents([''])
    }
  }, [bucket, reset])

  useEffect(() => {
    if (Object.keys(contentsSeen).length>0) { // Don't run at initialization
      setFirstRender(true)
      setContentsSeen({})
      setFetchContents([''])
      if (bucket == firstBucketPath[0]) {
        setPathArray(pathToArray(firstBucketPath[1], basePath));
      } else {
        setPathArray(pathToArray('', basePath));
      }
      if (reset) setReset(false); // Double-pass helps to fully clear
    }
  }, [reset])

  // Update main readout (path) whenever anything updates pathArray
  useEffect(() => {
    setPath(stripSlash(arrayToPath(pathArray, basePath)));
  }, [pathArray])

  // Detect when path has been updated externally
  useEffect(() => {
    const arrays_path = arrayToPath(pathArray, basePath);
    if (arrays_path != path && stripSlash(arrays_path) != path) {
      setPathArray(pathToArray(path, basePath));
    }
  }, [path]);

  // onChange for inners
  const updateTarget = (newPath: string, currentPathSet: string[], depth: number) => {
    setFirstRender(false)
    setShowAddButton(-1)
    
    // Also trim at the level just picked in case not actually the deepest
    let nextPathSet = [...currentPathSet].slice(0,depth+1)
    nextPathSet[depth] = newPath
    
    // Determine type of newPath target.  Can assume non-null answer because newPath has been chosen from known contents
    // Don't actually set here as calling useTargetType now likely triggers an infinite useEffect loop!
    const newType = (useDetermineType(
      newPath,
      stripSlash(arrayToPath(nextPathSet.slice(0,-1), basePath)),
      false
    ) as 'file' | 'folder')
    
    // Set new path value
    const targetPath = arrayToPath(nextPathSet, basePath)
    if (newType=='file') {
      setPathArray(nextPathSet)
    } else {
      // Also fetch folder contents (priority = here > fetchContents items) and add level if has usable contents
      if (!pathSeen(targetPath)) {
        fetchFolderContents(targetPath, (f) => {
          if (contentUse(f, allowFiles).length > 0) {
            nextPathSet = [...nextPathSet, '']
          }
          setPathArray(nextPathSet)
        })
      } else {
        if (contentUse(contentsSeen[pathInMetis(targetPath)] as folderSummary, allowFiles).length > 0) {
          nextPathSet = [...nextPathSet, '']
        }
        setPathArray(nextPathSet)
      }
    }

    useType(newType)
  }
  
  // clear for inners (a.k.a. select previous folder)
  const trimPath = (depth: number) => {
    const nextPathSet = depth > 0 ? [...pathArray].slice(0,depth) : ['']
    useType('folder')
    setPathArray(nextPathSet)
    if (depth > 0) setShowAddButton(depth-1)
  }

  return <Grid 
    key={`file-or-folder-picker-${project_name}-${bucket}`}
    container
    direction='column'
    className={className}
  >
    {pathArray.map( (p, index, fullArray) => {
      // const before = index > 1 ? fullSet[index-1] : '';
      const icon = index > 0 ? <SubdirectoryArrowRightIcon fontSize='small'/> : <i className="fas fa-folder" aria-hidden="true" style={{padding:'2px'}}/>;
      const thisPath = arrayToPath(fullArray.slice(0,index), basePath)
      const ready = pathSeen(thisPath)
      const thisLabel = (index==0) ? labelUse : undefined
      return (
        <Grid item container
          key={`sub-picker-${project_name}-${bucket}-`+ index}
          >
          <Grid item container alignItems='flex-end' style={{width: 'auto', paddingBottom: '7px', paddingRight: '2px'}}>
            {icon}
          </Grid>
          <Grid item style={{flex: '1 1 auto'}}>
            {ready ? <FileOrFolderInner
              optionSet={ contentUse(contentsSeen[pathInMetis(thisPath)] as folderSummary, allowFiles)}
              target={p}
              setTarget={(e: string) => updateTarget(e, pathArray, index)}
              onEmpty={ () => { trimPath(index) } }
              path={thisPath}
              autoOpen={!firstRender}
              placeholder={(index==0) ? topLevelPlaceholder : undefined}
              label={thisLabel}
              disablePortal={disablePortal}
            /> : <Autocomplete
              disabled
              options={['']}
              value=''
              disableClearable
              disablePortal={disablePortal}
              style={{paddingTop: thisLabel!==undefined ? 8 : undefined}}
              renderInput={(params) => (
                <TextField
                  {...params}
                  placeholder={'Awaiting contents or bucket choice'}
                  label={thisLabel}
                  size='small'
                  InputLabelProps={{shrink: true}}
                />
              )}
            />
            }
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
          {index!=0 ? null : <Grid item container alignItems='flex-end' style={{width: 'auto'}}>
            {/* <Tooltip title="Refresh Content Cache"> */}
              <IconButton
                aria-label="Refresh Content Cache"
                size="small"
                onClick={() => setReset(true)}>
                <AutorenewIcon/>
              </IconButton>
            {/* </Tooltip> */}
          </Grid>
          }
        </Grid>
      )
    })}
  </Grid>
}

export function PickFolder({ project_name=CONFIG.project_name, bucket, setPath, path, label, basePath, topLevelPlaceholder='', className, disablePortal=true}: {
  project_name?: string;
  bucket: string;
  setPath: Function;
  path: string; // the overall "value" / output
  label?: string;
  basePath: string; // an immutable portion of path.  Can be '' to access the entire bucket and to be compatible with exploring across buckets in sync with a PickBucket companion.
  topLevelPlaceholder?: string;
  className?: string;
  disablePortal?: boolean;
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
    className={className}
    disablePortal={disablePortal}
    />
}

