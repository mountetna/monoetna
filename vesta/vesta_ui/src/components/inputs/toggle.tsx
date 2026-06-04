import * as React from 'react';
import { styled } from '@mui/material';
import { Box } from '@mui/system';
import { useTheme } from '@mui/material';

const Toggle = ({active,setActive}:{
  active: boolean,
  setActive: (a:boolean) => void
}) => {
  return <Box
    onClick={ () => setActive(!active) }
    sx={{
      width: '49px',
      height: '27px',
      borderRadius: '27px',
      bgcolor: active ? 'green.grade100' : 'ground.grade75'
    }}>
    <Box
      sx={ theme => ({
        width: '21px',
        height: '21px',
        borderRadius: '21px',
        mt: '3px',
        ml: active ? '25px' : '3px',
        bgcolor: active ? 'ground.grade10' : 'ground.grade25',
        // @ts-ignore
        transition: theme.transitions.create(
          ['margin-left'],
          {
            // @ts-ignore
            easing: theme.transitions.easing.ease,
            // @ts-ignore
            duration: theme.transitions.duration.ease,
          }
        ),
      })}
    />
  </Box>;
}

export default Toggle;
