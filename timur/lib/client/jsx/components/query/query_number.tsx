import React, {useCallback, useMemo} from 'react';
import Grid from '@mui/material/Grid';
import Button from '@mui/material/Button';
import { makeStyles } from '@mui/styles';

const roman = (num) => {
  let roman = {
    m: 1000, cm: 900,
    d: 500, cd: 400,
    c: 100, xc: 90,
    l: 50, xl: 40,
    x: 10, ix: 9,
    v: 5, iv: 4, i: 1
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
    [ '1', false ],
    [ 'A', false ],
    [ 'i', true ],
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
