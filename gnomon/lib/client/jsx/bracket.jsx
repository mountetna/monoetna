import React from 'react';
import Grid from '@material-ui/core/Grid';

const Bracket = ({top, bottom, left, width}) => (
  <Grid
    style={{
      position: 'absolute',
      top, bottom, left, width,
      border: '1px solid black',
      height: '10px',
      [ top == undefined ? 'borderBottom' : 'borderTop' ]: 'none'
    }}
  />
)

export default Bracket;
