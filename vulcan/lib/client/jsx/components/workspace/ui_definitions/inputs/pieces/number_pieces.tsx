import React, { useEffect, useState } from 'react';
import TextField from '@material-ui/core/TextField';
import {IntegerInput, FloatInput} from 'etna-js/components/inputs/numeric_input.jsx'

function parseIntBetter(s: string) {
  const parsed = parseInt(s, 10);
  return parsed + '' == s ? parsed : NaN;
}

function from_num(num: number | null, float = false) {
  if (float) {
    return num == null ?
      '0.0' :
      num.toString().includes('.') ?
        num.toString() :
        num.toString() + '.0';
  }
  return num == null ? '' : num.toString();
}

// export function NumberPiece(
//   key: string, changeFxn: (value: number | null, key: string) => void, value: number | null = null,
//   label?: string,
//   minWidth: number = 200,
//   block_decimal: boolean = false
// ) {
//   const [inputState, setInputState] = useState({
//     text: from_num(value),
//     intVal: value,
//     hasError: false
//   });
//   // Detect external changes
//   useEffect(() => {
//     console.log("checking")
//     if (value != inputState.intVal) {
//       onNewNum({target: {value: value}});
//     }
//   }, [value, inputState.intVal])

//   function onNewNum(event: any) {
//     const parser = block_decimal ? parseIntBetter : Number;
//     const parsed = parser(event.target.value);
//     if (isNaN(parsed)) {
//       setInputState({text: event.target.value, intVal: null, hasError: true})
//       changeFxn(null, key);
//     } else {
//       setInputState({text: event.target.value, intVal: parsed, hasError: false});
//       changeFxn(parsed, key);
//     }
//   }
//   // console.log({inputState})

//   return (
//     <div style={{paddingTop: !!label ? 8 : 0}}>
//       <TextField
//         key={`${key}-number`}
//         value={inputState.text}
//         label={label}
//         error={inputState.hasError}
//         onChange={onNewNum}
//         size='small'
//         style={{minWidth: minWidth || 200}}
//       />
//     </div>
//   );
// }

export function NumberPiece(
  key: string, changeFxn: (value: number | null, key: string) => void, value: number | null = null,
  label?: string,
  minWidth: number = 200,
  integer: boolean = false
) {
  
  function onNewNum(event: any) {
    const parser = integer ? parseIntBetter : Number;
    const parsed = parser(event.target.value);
    if (!isNaN(parsed)) {
      changeFxn(parsed, key);
    } else {
      console.log(`Number update not used as it did not parse properly: ${event.target.value}`)
    }
  }
  
  return (
    <div style={{paddingTop: !!label ? 8 : 0}}>
      <TextField
        key={`${key}-number`}
        value={from_num(value, !integer)}
        label={label}
        onChange={onNewNum}
        size='small'
        style={{minWidth: minWidth || 200}}
      />
    </div>
  );
}

export function FloatPiece(
  key: string, changeFxn: (value: number | null, key: string) => void, value: number | null = null,
  label?: string,
  minWidth: number = 200
) {
  return NumberPiece(key, changeFxn, value, label, minWidth, false)
}

export function IntegerPiece(
  key: string, changeFxn: (value: number | null, key: string) => void, value: number | null = null,
  label?: string,
  minWidth: number = 200
) {
  return NumberPiece(key, changeFxn, value, label, minWidth, true)
}
