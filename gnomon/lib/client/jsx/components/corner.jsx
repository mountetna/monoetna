import React from 'react';
import Grid from '@material-ui/core/Grid';
import {makeStyles} from '@material-ui/core/styles';

const cornerStyle = makeStyles(theme => ({
  corner: {
    position: 'absolute',
    borderLeft: '1px solid black',
    opacity: 0.1,
    borderTop: '1px solid black'
  }
}));

const Corner = ({top, bottom, left, right, className}) => {
  const classes = cornerStyle();
  return <Grid
    className={ `${className} ${classes.corner}` }
    style={{
      bottom,
      left,
      height: top - bottom,
      width: right - left
    }}
  />;
};

export default Corner;
