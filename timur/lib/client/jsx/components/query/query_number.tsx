import React, {useCallback, useMemo} from 'react';
import Grid from '@mui/material/Grid';
import Button from '@mui/material/Button';
import { makeStyles } from '@mui/styles';

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

const useStyles = makeStyles((theme) => ({
  number: {
    minWidth: '25px',
    width: '25px',
    padding: '3.5px'
  }
}));

const QueryNumber = ({number, level, setRemoveHint, onClick}) => {
  const classes = useStyles();
  const numType = [
    [ 'I', true ],
    [ 'A', false ],
    [ '1', false ],
    [ 'a', false ]
  ];
  return <Button
    size='small'
    color='black'
    variant='variant'
    onMouseEnter={() => setRemoveHint && setRemoveHint(true)}
    onMouseLeave={() => setRemoveHint && setRemoveHint(false)}
    className={classes.number}
    style={{
      cursor: setRemoveHint ? 'pointer' : 'default'
    }}
    onClick={onClick}>
  {
    numType[level][1] ? roman(number+1) :
    String.fromCharCode(number + numType[level][0].charCodeAt(0))
  }.&nbsp;
</Button>
  }

export default QueryNumber;
