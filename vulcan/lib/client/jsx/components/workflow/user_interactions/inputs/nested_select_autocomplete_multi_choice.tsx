// Input component that takes nested object
//   and shows the keys one level at a time.
// Returns the last "Leaf" that the user selects.
import React, {useState, useEffect, useCallback, useMemo} from 'react';
import {WithInputParams} from './input_types';
import {some, Maybe} from '../../../../selectors/maybe';
import {useMemoized} from '../../../../selectors/workflow_selectors';
import {joinNesting} from './monoids';
import {useSetsDefault} from './useSetsDefault';
import { Grid } from '@material-ui/core'
import SubdirectoryArrowRightOutlinedIcon from '@material-ui/icons/SubdirectoryArrowRight';
import { arrayLevels } from './user_input_pieces';
import SelectAutocompleteMultiPickInput from './select_autocomplete_multi_choice';

type OptionSet = {[k: string]: null | OptionSet};

export const sep = '---'

export function flattenOptionPaths(options: OptionSet, pre = [] as string[]): {[key:string]: string[]} {
  // Output: keys = 'value' options, values = array holding path of keys upstream in OptionSet leading to said 'value' 
  let paths = {} as {[key:string]: string[]}
  for (let [key, value] of Object.entries(options)) {
    if (value == null) {
      paths[key] = pre
    } else {
      paths = {...paths, ...flattenOptionPaths(value, pre.concat([key]))}
    }
  }
  return paths
}

export function leafParentPaths(pathMap: ReturnType<typeof flattenOptionPaths>, _sep: string = sep): string[] {
  // Output: All "Option Paths", unique paths (now combined from [level1, level2].join(sep) which hold leaf options
  return arrayLevels(Object.values(pathMap).map(val => val.join(_sep)))
}

function targettedPathValues(values: string[] | null, pathMap: ReturnType<typeof flattenOptionPaths>): {[key:string]: string[]} {
  // Output: keys = Used "Option Paths", a.k.a. ones that a 'value' was selected from; values = selected 'value's from those paths.
  if (values == null || values.length<1 || Object.keys(pathMap).length<1) return {}
  let pathVals = {} as {[key:string]: string[]}
  for (let value of values) {
    const vals_path = pathMap[value].join(sep)
    if (Object.keys(pathVals).includes(vals_path)) {
      pathVals[vals_path] = [...pathVals[vals_path], value]
    } else {
      pathVals[vals_path] = [value]
    }
  }
  return pathVals
}

export function pathValues(pathString: string, allOptions: OptionSet, _sep: string = sep): string[] {
  // Output: All 'value' options which come from this path
  const pathLevels = pathString.split(_sep)
  let values = {...allOptions}
  pathLevels.forEach( (this_level) => {
    values = values[this_level] as OptionSet
  })
  return Object.keys(values)
}

function valuesPerSelectPaths(valuesPerPath: {[key:string]: string[]}, paths: string | string[]) {
  // Output: selected 'value's which come from this / these path(s)
  if (!Array.isArray(paths)) {
    paths = [paths]
  }
  const output = {} as typeof valuesPerPath
  for (let path of paths) {
    output[path] = valuesPerPath[path]
  }
  return output
}

function valuesFromValuesPerPath(valuesPerPath: {[key:string]: string[]}) {
  // Output: single array containing all chosen 'value's / all values from a 'targettedPathValues()' or 'valuesPerSelectPaths()' output
  return arrayLevels(([] as string[]).concat.apply([] as string[], Object.values(valuesPerPath))) as string[]
}

export default function NestedSelectAutocompleteMultiPickInput({ label, testIdAppend='', data, onChange, ...props }: WithInputParams<{
  label?: string
  testIdAppend?: string
}, string[], OptionSet>) {
  const value: string[] = useSetsDefault([] as string[], props.value, onChange);
  const allOptions = useMemoized(joinNesting, data);
  const [valuesPerPath, setValuesPerPath] = useState({} as {[key:string]: string[]});
  const [paths, setPaths] = useState([] as string[]);

  const flattenedPaths = useMemo(() => {
    return flattenOptionPaths(allOptions)
  }, [allOptions])

  const pathOptions = useMemo(() => {
    return leafParentPaths(flattenedPaths)
  }, [flattenedPaths])

  useEffect(() => {
    if (value!=null) {
      const vals_per_path = targettedPathValues(value, flattenedPaths)
      const needed_paths = arrayLevels(Object.keys(vals_per_path))
      setValuesPerPath(vals_per_path)
      // Re-add missing paths (likely from page reload)
      if ( needed_paths.length>0 && ! needed_paths.every( (path) => paths.includes(path) )) {
        setPaths(arrayLevels(paths.concat(Object.keys(vals_per_path))));
      }
    }
  }, [flattenedPaths, value]);

  // OptionPath Selection
  const handlePickedOptionPaths = useCallback( (event: any, e: string[]) => {
    const e_use = arrayLevels(e)
    const currentPaths = [...paths]
    if (e_use != currentPaths) {
      const prevPaths = paths.filter((path:string) => e_use.includes(path))
      setPaths(e_use)
      // Remove associated values if user removed an option path
      // (prevPaths = chosen paths that had previously been chosen)
      if (currentPaths.length > prevPaths.length) {
        onChange(some(valuesFromValuesPerPath(valuesPerSelectPaths(valuesPerPath, prevPaths))))
      }
    }
  }, [valuesPerPath])
  const optionPathPicker = <Grid item>
      <SelectAutocompleteMultiPickInput
        key='MultiMultiPick-optionPaths'
        onChange={(v) => {}}
        onChangeOverride={handlePickedOptionPaths}
        value={some(paths)}
        data={{a: pathOptions}}
        placeholder={paths.length<1 ? 'Option Sets' : undefined}
        label={label}
        testId={testIdAppend ? `MultiMultiPick-optionPaths-${testIdAppend}`: 'MultiMultiPick-optionPaths'}
      />
    </Grid>

  // Value Selection
  const handlePickedValues = useCallback( (event: any, e: string[], pathString: string) => {
    const e_use = arrayLevels(e)
    const newValuesPerPath = {...valuesPerPath}
    newValuesPerPath[pathString] = e_use
    const newValues = valuesFromValuesPerPath(newValuesPerPath)
    onChange(some(newValues))
  }, [valuesPerPath])
  const singleValuePicker = (pathString: string, values: string[], disabled = false) => {
    return(
      <Grid item container direction='row' key={pathString} alignItems='flex-start'>
        <Grid item xs={1} container justifyContent='flex-end'>
          <SubdirectoryArrowRightOutlinedIcon fontSize='small' color="secondary" style={{paddingTop: '3px', paddingRight: '2px'}}/>
        </Grid>
        <Grid item xs={11}>
          <SelectAutocompleteMultiPickInput
            key={pathString+'-leaves'}
            data={disabled ? {a: ['']} : {a: pathValues(pathString, allOptions)}}
            disabled={disabled}
            placeholder={disabled ? 'Awaiting Option Set selection' : undefined}
            label={pathString+' options'}
            value={some(values)}
            onChangeOverride={ (event: any, e: string[]) => handlePickedValues(event, e, pathString) }
            onChange={(v) => {}}
            testId={testIdAppend ? `MultiMultiPick-${pathString}-leaves-${testIdAppend}`: `MultiMultiPick-${pathString}-leaves`}
          />
        </Grid>
      </Grid>
    )}
  const allValuePickers = paths.length==0 ? singleValuePicker('TBD', [], true) : <Grid item>
      {paths.map( (pathString) => {
        const values_use = Object.keys(valuesPerPath).includes(pathString) ? valuesPerPath[pathString] : []
        return singleValuePicker(pathString, values_use)
      } )}
    </Grid>

  return <Grid container direction='column'>
    {optionPathPicker}
    {allValuePickers}
  </Grid>
};
