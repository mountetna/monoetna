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

const sep = '---'

function flattenOptionPaths(options: OptionSet, pre = [] as string[]): {[key:string]: string[]} {
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

function leafParentPaths(pathMap: ReturnType<typeof flattenOptionPaths>, _sep: string = sep): string[] {
  // Output: All "Option Paths", unique paths (now combined from [level1, level2].join(sep) which hold leaf options
  return arrayLevels(Object.values(pathMap).map(val => val.join(_sep)))
}

function targettedPathValues(values: string[] | null, pathMap: ReturnType<typeof flattenOptionPaths>): {[key:string]: string[]} {
  // Output: keys = Used "Option Paths", a.k.a. ones that a 'value' was selected from; values = selected 'value's from those paths.
  if (values == null || Object.keys(pathMap).length<1) return {}
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

function pathValues(pathString: string, allOptions: OptionSet, _sep: string = sep): string[] {
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
  return arrayLevels([].concat.apply([], Object.values(valuesPerPath))) as string[]
}

export default function NestedSelectAutocompleteMultiPickInput({ label, data, onChange, ...props }: WithInputParams<{
  label?: string
}, string[], OptionSet>) {
  const value = useSetsDefault(null, props.value, onChange) as Maybe<string[]>;
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
    if (value!=null && value.length>0) {
      const vals_per_path = targettedPathValues(value, flattenedPaths)
      const needed_paths = arrayLevels(Object.keys(vals_per_path))
      setValuesPerPath(vals_per_path)
      // Re-add missing paths (likely from page reload)
      if ( ! needed_paths.every( (path) => paths.includes(path) )) {
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
      // (prevPaths = chosen paths that would exist in valuesPerPath)
      if (prevPaths.length>0 && currentPaths != prevPaths) {
        onChange(some(valuesFromValuesPerPath(valuesPerSelectPaths(valuesPerPath, prevPaths))))
      }
    }
  }, [valuesPerPath])
  const optionPathPicker = <Grid item>
      <SelectAutocompleteMultiPickInput
        key={'MultiMultiPick-optionPsths'}
        onChange={(v) => {}}
        onChangeOverride={handlePickedOptionPaths}
        value={some(paths)}
        data={{a: pathOptions}}
        label={label? label+' - Option Paths' : 'Option Paths'}
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
  const singleValuePicker = (pathString: string, values: string[], disabled = false) => <Grid item container direction='row' key={pathString}>
    <Grid item container alignItems='flex-start' style={{width: 'auto', paddingTop: '3px', paddingRight: '2px'}}>
      <SubdirectoryArrowRightOutlinedIcon fontSize='small' color="secondary"/>
    </Grid>
    <Grid item style={{flex: '1 1 auto'}}>
      <SelectAutocompleteMultiPickInput
        key={pathString+'-leaves'}
        data={disabled ? {a: ['']} : {a: pathValues(pathString, allOptions)}}
        disabled={disabled}
        placeholder={disabled ? 'Awaiting Option Set selection' : undefined}
        label={pathString+' choices'}
        value={some(values)}
        onChangeOverride={ (event: any, e: string[]) => handlePickedValues(event, e, pathString) }
        onChange={(v) => {}}
      />
    </Grid>
  </Grid>
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
