'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Image from '../image/image';

const ThemeIcon = ({ theme }:{
  theme: Theme
}) => {
  return <Box
      sx={{
          display: 'flex',
          p: '6px',
          bgcolor: theme.color,
          borderRadius: '50%',
      }}
  >
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
}

export default ThemeIcon;
