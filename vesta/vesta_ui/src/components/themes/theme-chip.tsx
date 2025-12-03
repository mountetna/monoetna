'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Image from '../image/image';
import Typography from '@mui/material/Typography';

const ThemeChip = ({ theme }:{
  theme: Theme
}) => {
  return <Box
      sx={{
          display: 'flex',
          p: '6px',
          alignItems: 'center',
          bgcolor: theme.color,
          borderRadius: '25px',
      }}
  >
    <Box sx={{ display: 'flex', alignItems: 'center', px: '3px' }}>
      <Image
          src={theme.icon}
          alt={`Abstract icon for "${theme.name}" theme`}
          width={20}
          height={20}
          style={{
              backgroundColor: theme.color,
          }}
          hideBeforeLoad={false}
      />
    </Box>
    <Box sx={{ display: 'flex', alignItems: 'center', px: '3px' }}>
      <Typography sx={{ color: 'ground.grade10', textTransform: 'uppercase', whiteSpace: 'nowrap' }} variant="pBodyMono">
        { theme.name }
      </Typography>
    </Box>
  </Box>
}

export default ThemeChip;
