import React, { useEffect } from 'react';
import InputLabel from '@material-ui/core/InputLabel';
import Slider from '@material-ui/core/Slider';
import Grid from '@material-ui/core/Grid';
import { FloatPiece } from './number_pieces';
import { DataEnvelope } from '../../input_types';
import DropdownPiece, { DropdownPieceRct } from './dropdown_piece';
import TextField from '@material-ui/core/TextField';

export function val_wrap(v: any): DataEnvelope<typeof v> {
  return {'a': v};
}

export function key_wrap(k: string[]) {
  let de: DataEnvelope<null> = {};
  for (let ind = 0; ind < k.length; ind++) {
    de[k[ind]]=null;
  }
  return de;
}

export function arrayLevels(original: any[]) {
  function onlyUnique(value: any, index: number, self: any) {
    return self.indexOf(value) === index;
  }
  return Array.from(original).filter(onlyUnique);
}

export type PieceBaseInputs<V> = {
  name: string;
  changeFxn: (v: V, k?: string) => void;
  value: V;
  label?: string;
};

export function DisabledTextbox({
  name,
  label,
  text
}: {
  name: string
  label: string | undefined
  text: string
}): React.ReactElement {
  return <div key={name}>
    <TextField
      key={name}
      label={label}
      InputLabelProps={{ shrink: true }}
      size="small"
      disabled
      fullWidth
      placeholder={text}
    />
  </div>
};

/*
"Pieces" which follow define components/elements which can be used to fill in discrete parts of a user input widget.
They are named based on types of input methods.
Overall Output Structure Requirement:
  - It is assumed that the overarching widget will produce an output consisting of a hash of key/value pairs.
Some "pieces" have additional inputs, but the first 4 are ALWAYS:
  - key = the output key or name of the value element that this component should fill in. Also used as the component's key for the DOM
  - changeFxn = a function, which should be defined inside the larger ui-component in which these pieces are used.
      This funciton should take in 1) the new value and 2) the target 'key' / value elements' name, then utilize onChange to log the update.
  - value = the current value of this target element
  - label = a text label to be displayed with this component
*/

export function ReductionSetupPiece(
  key: string, changeFxn: Function, value: (string|null)[] = [null, '1', '2'],
  label: string[] = ['Dimensionality Reduction (DR)', 'x-axis DR Compenent', 'y-axis DR Component'],
  reduction_opts: DataEnvelope<string[]>) {
  
  function changeReduction(newElement: string|null) {
    return [newElement, '1', '2']
  }
  function changeDim(newElement: string|null, dim: 1|2) {
    let newValue = [...value]
    newValue[dim] = newElement
    return newValue
  }

  useEffect(()=>{
    // Default to _Recommended_ (which comes from 'umap_name' attribute of the dataset record)
    if (value[0]==null && reduction_opts != null && Object.keys(reduction_opts).includes('_Recommended_')) {
      changeReduction('_Recommended_')
    }
  }, [value, reduction_opts])

  const error_with_input = !reduction_opts;

  const disable_dims = value[0]==null || error_with_input;
  // console.log({reduction_opts})
  return(
    <Grid 
      key={key}
      container
      direction='column'
    >
      <Grid item>
        <DropdownPieceRct
          name={key+'-reduction'}
          changeFxn={(newElement: string | null) => changeFxn(changeReduction(newElement), key)}
          value={error_with_input ? 'Error in input setup' : value[0]}
          label={label[0]}
          options_in={error_with_input ? ['Error in input setup'] : Object.keys(reduction_opts)}
          minWidth={200}
          disabled={error_with_input}
        />
      </Grid>
      {[{i:1,ax:'x'},{i:2,ax:'y'}].map((v) => {
        const i: number = v.i;
        const k = `${key}-dim${v.ax}`
        return <Grid item key={k} style={{paddingLeft: 10}}>
          <DropdownPieceRct
            name={k}
            changeFxn={(newElement: string | null) => changeFxn(changeDim(newElement, i as 1 | 2), key)}
            value={value[i]}
            label={label[i]}
            options_in={value[0]==null ? ['1', '2'] : reduction_opts[value[0]] || [`${i}`]}
            minWidth={200}
            disabled={disable_dims}
          />
        </Grid>
      })}
    </Grid>
  )

}
