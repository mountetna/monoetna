import React, {Component} from 'react';
import Box from '@material-ui/core/Box';

const Legend = ({labels, width, height, onClick}) =>
  <Box sx={{
    display: 'flex',
    justifyContent: 'center',
    width, height,
    '& .legend_box': {
      display: 'inline-block',
      marginRight: '8px',
      width: '12px',
      height: '10px',
    },
    '& .legend_text': {
      display: 'inline',
      fontSize: '14px'
    },
    '& .legend_item': {
      cursor: 'pointer',
      display: 'flex',
      height: '100%',
      marginRight: '24px',
      alignItems: 'center',
      '&:hover': {
        '& .legend_box': {
          border: '2px solid black'
        }
      },
    }
  }}>
    {
      labels.flat().map(({name,color,disabled}) =>
        <Box key={ name } className='legend_item' onClick={ () => onClick(name)} >
          <Box className='legend_box' sx={{ background: disabled ? 'none' : color, border: `2px solid ${color}` }} />
          <Box className='legend_text'>{name}</Box>
        </Box>
      )
    }
  </Box>;

export default Legend;
