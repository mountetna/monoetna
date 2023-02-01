import React, { useState, useCallback, useEffect, useMemo } from 'react';
import {metisPath} from 'etna-js/api/metis_api.js'
import Autocomplete from '@material-ui/lab/Autocomplete';
import { TextField } from '@material-ui/core';
import { json_get } from 'etna-js/utils/fetch';
import SubdirectoryArrowRightIcon from '@material-ui/icons/SubdirectoryArrowRight';

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

function simplifyFolderListReturn(fullList: any) {
  return {
    folders: fullList.folders.map( (f: any) => f.folder_name) as string[],
    files: fullList.files.map( (f: any) => f.file_name) as string[]
  }
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

function SingleFileOrFolder({ project_name, bucket, path, setTarget, onEmpty, allowFiles = true}: {
  project_name?: string;
  bucket: string;
  path: string;
  setTarget: Function;
  onEmpty: Function;
  allowFiles?: boolean;
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
                folder_metis_return => {
                  if (Object.values(folder_metis_return).flat().length > 0 ) {
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
      options={targetList.concat('')}
      autoHighlight
      value={target}
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
      disableClearable={false}
      disablePortal
      inputValue={inputState}
      onInputChange={(event: any, newInputState: string) => {
        setInputState(newInputState);
      }}
      // style={{minWidth: 100}}
      renderInput={(params) => (
        <TextField
          {...params}
          error={ inputState != nullToEmptyString(target) }
          autoFocus
          size='small'
          InputLabelProps={{shrink: true}}
        />
      )}
    />
  );
}

export function PickFileOrFolder({ project_name=CONFIG.project_name, bucket, setPath, path, allowFiles}: {
  project_name?: string;
  bucket: string;
  setPath: Function;
  path: string;
  allowFiles?: boolean
}){
  const [pathArray, setPathArray] = useState(pathToArray(nullToEmptyString(path)));

  useEffect(() => {
    setPathArray(pathToArray(nullToEmptyString(path)));
  }, [path])

  // onChange for inners
  const updateTarget = useCallback(
    (newPath: string, currentPathSet: string[], depth: number) => {
      // Also trim at the level just picked in case not actually the deepest
      // console.log(newPath)
      let nextPathSet = [...currentPathSet].slice(0,depth+1)
      nextPathSet[depth] = newPath
      setPath(arrayToPath(nextPathSet))
    }, [setPath])
  
  // clear for inners
  const trimPath = useCallback(
    (depth: number) => {
    const nextPathSet = depth > 0 ? [...pathArray].slice(0,depth+1) : ['']
    setPath(arrayToPath(nextPathSet))
  }, [setPath])

  const targetSelectors = useMemo(() => {
    return pathArray.map( (p, index, fullSet) => {
      const before = index > 1 ? pathArray[index-1] : '';
      const indent = index > 1 ? (index-1)*10 : 0;
      const icon = index > 0 ? <SubdirectoryArrowRightIcon/> : null;
      return (
        <div 
          style={{display: 'inline-flex', paddingLeft: indent}}
          key={`sub-picker-${bucket}-${project_name}-${before}`+ index}
          >
          {icon}
          <SingleFileOrFolder
            project_name={project_name}
            bucket={bucket}
            setTarget={(e: string) => updateTarget(e, pathArray, index)}
            onEmpty={ () => { trimPath(index) } }
            path={arrayToPath(fullSet.slice(0,index))}
            allowFiles={allowFiles}
          />
        </div>
      )
    })
  }, [pathArray, project_name, bucket, updateTarget])

  // console.log(path)
  // console.log(pathArray)

  return bucket == null ? null : <div key={`file-or-folder-picker-${bucket}-${project_name}`}>
    {targetSelectors}
  </div>
}


