import React from 'react';
import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
import {makeStyles} from '@material-ui/core/styles';

const useStyles = makeStyles((theme) => ({
  letter: {
    fontSize: '30px'
  },
  box: {
    position: 'relative',
    width: '40px',
    height: '40px'
  },
}));

const Letter = ({letter, className}) => {
  const classes = useStyles();
  return <Grid className={`${classes.box} ${className}`}
    container alignItems='center' justify='center'>
    <span className={classes.letter}>{letter}</span>
  </Grid>
}

export default Letter;
