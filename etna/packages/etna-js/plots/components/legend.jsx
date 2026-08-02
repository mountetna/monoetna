import React, {Component} from 'react';
import Box from '@material-ui/core/Box';

const Legend = ({labels, width, height}) =>
  <Box sx={{
    display: 'flex',
    justifyContent: 'center',
    width, height
  }}>
    {
      labels.flat().map(({name,color}) =>
        <Box key={ name } sx={{
          display: 'flex',
          height: '100%',
          marginRight: '24px',
          alignItems: 'center',
        }}>
          <Box sx={{
            background: color,
            display: 'inline-block',
            marginRight: '8px',
            width: '16px',
            height: '8px',
            verticalAlign: 'middle',
          }}/>
          <Box sx={{
            display: 'inline',
            fontSize: '12px'
          }}>{name}</Box>
        </Box>
      )
    }
  </Box>;

export default Legend;
