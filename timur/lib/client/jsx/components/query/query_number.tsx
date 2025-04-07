import React, {useCallback, useMemo} from 'react';
import Grid from '@material-ui/core/Grid';
import Button from '@material-ui/core/Button';
import { makeStyles } from '@material-ui/core/styles';

const roman = (num) => {
  let roman = {
    m: 1000, cm: 900,
    d: 500, cd: 400,
    c: 100, xc: 90,
    l: 50, xl: 40,
    x: 10, ix: 9,
    v: 5, iv: 4, i: 1
  };

  num = num + 1;

  return Object.keys(roman).reduce( (str, i) => {
    let q = Math.floor(num / roman[i]);
    num -= q * roman[i];
    return str + i.repeat(q);
  }, '');
}

const letter = (num) => {
  return String.fromCharCode((num % 26) + 'A'.charCodeAt(0)).repeat( Math.floor(num / 26) + 1 );
}

const numeric = (num) => {
  return String(num + 1);
}

const useStyles = makeStyles((theme) => ({
  number: {
    minWidth: '25px',
    width: '25px',
    padding: '3.5px'
  }
}));

const numType = [
  { type: numeric },
  { type: letter, caps: true },
  { type: roman, caps: false },
  { type: letter, caps: false }
];

const QueryNumber = ({number, level, setRemoveHint, onClick}) => {
  const classes = useStyles();

  let formattedNumber = numType[level].type(number);

  formattedNumber = numType[level].caps ? formattedNumber.toUpperCase() : formattedNumber;

  return <Button
    size='small'
    variant='text'
    onMouseEnter={() => setRemoveHint && setRemoveHint(true)}
    onMouseLeave={() => setRemoveHint && setRemoveHint(false)}
    className={classes.number}
    color={
      setRemoveHint ? 'primary' : 'default'
    }
    style={{
      cursor: setRemoveHint ? 'pointer' : 'text'
    }}
    onClick={onClick}>
    { formattedNumber }.&nbsp;
  </Button>
}

export default QueryNumber;
