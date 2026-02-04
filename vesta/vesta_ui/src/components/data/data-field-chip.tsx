'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';

const DataFieldChip = ({ value, count }:{
  value: string,
  count?: number
}) => {
  return <Box
      sx={{
          display: 'flex',
          py: '4px',
          px: '8px',
          alignItems: 'center',
          bgcolor: 'ground.grade75',
          borderRadius: '8px',
      }}
  >
    { count !== undefined && 
      <Typography sx={{
        color: 'ground.grade10',
        bgcolor: 'ground.grade100',
        p: '6px',
        mr: '10px',
        borderRadius: '50px'
      }} variant="pBodyMono">
        { count }
      </Typography>
    }
      <Typography sx={{ color: 'ground.grade10' }} variant="pBodyMono">
        { value }
      </Typography>
  </Box>
}

export default DataFieldChip;
