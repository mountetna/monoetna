import React from 'react';
import {makeStyles} from '@material-ui/core/styles';
import Grid from '@material-ui/core/Grid';

const lineStyles = makeStyles((theme) => ({
  line: {
    transformOrigin: 'bottom left',
    height: 0
  }
}));

const Line = ({bottom, left, top, right, className}) => {
  const classes = lineStyles();

  const h = top - bottom;
  const w = left - right;
  const r = Math.sqrt(h*h + w*w);
  const a = Math.atan(h/w);

  return <Grid
    className={ `${className} ${classes.line}` }
      style={{
        left, bottom,
        width: r,
        transform: `rotate(${a}rad)`,
    }}
  />;
};

export default Line;
