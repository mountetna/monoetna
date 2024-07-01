import React, {useCallback, useMemo} from 'react';
import Grid from '@mui/material/Grid';
import Typography from '@mui/material/Typography';

const roman = (num) => {
  let roman = {
    M: 1000, CM: 900,
    D: 500, CD: 400,
    C: 100, XC: 90,
    L: 50, XL: 40,
    X: 10, IX: 9,
    V: 5, IV: 4, I: 1
  };

  return Object.keys(roman).reduce( (str, i) => {
    let q = Math.floor(num / roman[i]);
    num -= q * roman[i];
    return str + i.repeat(q);
  }, '');
}

const QueryNumber = ({number, level}) => {
  const numType = [
    [ 'I', true ],
    [ 'A', false ],
    [ '1', false ],
    [ 'a', false ]
  ];
  return <Typography>
  {
    numType[level][1] ? roman(number+1) :
    String.fromCharCode(number + numType[level][0].charCodeAt(0))
  }.&nbsp;
</Typography>
  }

export default QueryNumber;
