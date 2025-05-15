// Input component that takes nested object
//   and shows the keys one level at a time.
// Returns the last "Leaf" that the user selects.
import React, {useState, useEffect, useCallback, useMemo} from 'react';
import { Grid } from '@material-ui/core'
import SubdirectoryArrowRightOutlinedIcon from '@material-ui/icons/SubdirectoryArrowRight';
import { arrayLevels, key_wrap, PieceBaseInputs } from './user_input_pieces';
import { DropdownMultiChoicePieceRct } from './dropdown_multi_choice_piece';
import { nestedOptionSet, OptionSet } from '../../input_types';

export const sep = '---'

export function flattenOptionPaths(options: nestedOptionSet, pre = [] as string[]): {[key:string]: string[]} {
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

export function pathValues(pathString: string, allOptions: nestedOptionSet, _sep: string = sep): string[] {
  // Output: All 'value' options which come from this path
  const pathLevels = pathString.split(_sep)
  let values = {...allOptions}
  pathLevels.forEach( (this_level) => {
    values = values[this_level] as nestedOptionSet
  })
  return Object.keys(values).filter(k => values[k]==null)
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

export interface NestedDropdownMultiChoicePieceInputs extends PieceBaseInputs<string[] | null> {
  options_in: OptionSet,
  testId?: string
  sorted?: boolean
}

export function NestedDropdownMultiChoicePiece(
  key: NestedDropdownMultiChoicePieceInputs['name'],
  changeFxn: NestedDropdownMultiChoicePieceInputs['changeFxn'],
  value: NestedDropdownMultiChoicePieceInputs['value'],
  label: NestedDropdownMultiChoicePieceInputs['label'] = '',
  options: string[] | nestedOptionSet | null,
  sorted: boolean = true
): React.ReactElement | null {
  if (options==null) return null;

  return <NestedDropdownMultiChoicePieceRct
    name={key}
    label={label}
    testId={key}
    value={value}
    options_in={options}
    changeFxn={changeFxn}
  />
}

export default function NestedDropdownMultiChoicePieceRct({
  name,
  changeFxn,
  value,
  label,
  options_in,
  testId = undefined,
  sorted = false
}: NestedDropdownMultiChoicePieceInputs) {
  const [valuesPerPath, setValuesPerPath] = useState({} as {[key:string]: string[]});
  const [paths, setPaths] = useState([] as string[]);
  const options: nestedOptionSet = useMemo(() => {
    if (Array.isArray(options_in)) {
      return key_wrap([...options_in]) as nestedOptionSet
    }
    return options_in
  }, [options_in])

  const flattenedPaths = useMemo(() => {
    return flattenOptionPaths(options)
  }, [options])

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
  const handlePickedOptionPaths = useCallback( (v: string[] | null, k?:string) => {
    const v_use = arrayLevels(v || [])
    const currentPaths = [...paths]
    if (v_use != currentPaths) {
      const prevPaths = paths.filter((path:string) => v_use.includes(path))
      setPaths(v_use)
      // Remove associated values if user removed an option path
      // (prevPaths = chosen paths that had previously been chosen)
      if (currentPaths.length > prevPaths.length) {
        changeFxn(valuesFromValuesPerPath(valuesPerSelectPaths(valuesPerPath, prevPaths)), k)
      }
    }
  }, [valuesPerPath])
  const optionPathPicker = <Grid item>
      <DropdownMultiChoicePieceRct
        name='MultiMultiChoice-optionPaths'
        changeFxn={(v: string[] | null, k?: string) => handlePickedOptionPaths(v, name)}
        value={paths}
        options_in={pathOptions}
        placeholder={paths.length<1 ? 'Option Sets' : undefined}
        label={label}
        testId={!!testId ? `MultiMultiChoice-optionPaths-${testId}`: undefined}
      />
    </Grid>

  // Value Selection
  const handlePickedValues = useCallback( (pathString: string, v: string[] | null, k?: string) => {
    const v_use = arrayLevels(v || [])
    const newValuesPerPath = {...valuesPerPath}
    newValuesPerPath[pathString] = v_use
    const newValues = valuesFromValuesPerPath(newValuesPerPath)
    changeFxn(newValues, k)
  }, [valuesPerPath])
  const singleValuePicker = (pathString: string, values: string[], disabled = false) => {
    return(
      <Grid item container direction='row' key={pathString} alignItems='flex-start'>
        <Grid item xs={1} container justifyContent='flex-end'>
          <SubdirectoryArrowRightOutlinedIcon fontSize='small' color="secondary" style={{paddingTop: '3px', paddingRight: '2px'}}/>
        </Grid>
        <Grid item xs={11}>
          <DropdownMultiChoicePieceRct
            name={pathString+'-leaves'}
            options_in={disabled ? [''] : pathValues(pathString, options)}
            disabled={disabled}
            placeholder={disabled ? 'Awaiting Option Set selection' : undefined}
            label={pathString+' options'}
            value={values}
            changeFxn={ (v: string[] | null, k?: string) => handlePickedValues(pathString, v, name) }
            testId={!!testId ? `MultiMultiChoice-${pathString}-leaves-${testId}`: undefined}
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

  return <Grid key={name} container direction='column'>
    {optionPathPicker}
    {allValuePickers}
  </Grid>
};
