import React, {useState, useEffect, useCallback} from 'react';
import DropdownPiece from './dropdown_piece';
import { key_wrap } from './user_input_pieces';

type NestedOptionSet = {[k: string]: null | NestedOptionSet};

function getPath(options: NestedOptionSet, leaf: string): string[] {
  for (let [key, value] of Object.entries(options)) {
    if (value && typeof value === 'object') {
      // Look one step ahead for our leaf nodes
      if (leaf in value && null == value[leaf])
        {return [key, leaf];}
      const path = getPath(value, leaf);
      if (path.length > 0) return [key, ...getPath(value, leaf)];
    } else if (null == value && key == leaf) {
      return [key];
    }
  }

  return [];
}

function getOptions(
  desiredPath: string[],
  optionSet: NestedOptionSet
): string[] | null {
  if (null == desiredPath || desiredPath.length === 0)
    {return Object.keys(optionSet);}
  // _.at always returns an array, so unwrap it.
  let results = _.at(optionSet, desiredPath.join('.'))[0];
  if (results) return Object.keys(results);
  return null;
}

type NestedDropdownPieceInputs = {
  name?: string,
  changeFxn: (v: string | null, k?: string) => void, // k likely becomes meaningless so should not be used in the function
  value: string | null,
  label?: string,
  nestedOptions: NestedOptionSet,
  disabled?: boolean
}

export default function NestedDropdownPiece(
  key: NestedDropdownPieceInputs['name'],
  changeFxn: NestedDropdownPieceInputs['changeFxn'],
  value: NestedDropdownPieceInputs['value'],
  label: NestedDropdownPieceInputs['label'] = '',
  options: NestedDropdownPieceInputs['nestedOptions'] | string[],
  sorted: boolean = true,
  minWidth: number = 200,
  disabled: NestedDropdownPieceInputs['disabled'] = false,
): React.ReactElement | null {
  if (options==null) return null;
  if (Array.isArray(options)) {
    options = key_wrap([...options]) as {[k:string]: null}
  }
  return <NestedDropdownPieceRct
    name={key}
    changeFxn={changeFxn}
    value={value}
    label={label}
    nestedOptions={options}
    disabled={disabled}
  />
}

export function NestedDropdownPieceRct({
  name,
  changeFxn,
  value,
  label = '',
  nestedOptions,
  disabled = false,
}: NestedDropdownPieceInputs): React.ReactElement | null {
  const [path, setPath] = useState([] as string[]);
  const [leafOptions, setLeafOptions] = useState(Object.keys(nestedOptions) as string[] | null);

  useEffect(() => {
    if (value!=null && path.length==0) {
      const updatedPath = getPath(nestedOptions, value);
      setPath(updatedPath);
      setLeafOptions(getOptions(updatedPath, nestedOptions));
    }
  }, [nestedOptions, value]);

  const handleSelect = useCallback(
    (value: string | null, depth: number) => {
      // User has not selected something...perhaps
      //   still typing?
      if (null == value) return;

      // figure out where the next selection comes from, check if we're picking a leaf node.
      const updatedPath = path.slice(0, depth);
      updatedPath.push(value);
      const updatedLeafOptions = getOptions(updatedPath, nestedOptions);
      setLeafOptions(updatedLeafOptions);
      setPath(updatedPath);

      if (updatedLeafOptions == null) {
        // If we are updating a leaf
        changeFxn(value, name);
      } else {
        // Otherwise a leaf has not been selected.
        changeFxn(null, name);
      }
    },
    [nestedOptions, path, changeFxn]
  );

  // console.log({value})
  // console.log({path})
  return (
    <div key='nested-dropdown-input'>
      {path.map((value, index) => {
        const options = getOptions(path.slice(0, index), nestedOptions);
        return (
          !!options && // skip here if leaf
          <div key={`nested-dropdown-${index}`}>
            {DropdownPiece(
              `nested-dropdown-${index}`,
              (v, k) => handleSelect(v, index),
              value,
              index==0 ? label : undefined,
              options,
              false,
              200,
              disabled
            )}
          </div>
        );
      })}
      { !!leafOptions &&
        (value==null || leafOptions.includes(value)) &&
        DropdownPiece(
          `nested-dropdown-${path.length}`,
          (v, k) => handleSelect(v, path.length),
          value,
          path.length == 0 ? label : undefined,
          leafOptions,
          false,
          200,
          disabled
        )
      }
    </div>
  );
};
