import React from 'react';
import Grid from '@material-ui/core/Grid';

const Bracket = ({top, bottom, left, width, className}) => (
  <Grid
    className={className}
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
