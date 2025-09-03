'use client'

import * as React from 'react';
import Box from '@mui/system/Box';
import Grid from '@mui/system/Grid';
import Button from '@mui/material/Button';
import Typography from '@mui/material/Typography';
import Image from '../../components/image/image';
import projectExplorerBackshape from '/public/images/project-explorer-backshape.svg'

export default function ProjectInfo() {
  return <Box sx={{ padding: '90px 144px' }}>
      <Typography
          variant='h2'
      >
          Projects
      </Typography>
      <Box sx={{ position: 'relative' }}>
        <Grid container spacing={8}>
          <Typography variant="pLarge">
              The Data Library hosts a series of research projects containing multi-omic data and de-identified clinical data. Project data is shared in our research community to best leverage these rich datasets for biological discovery.
          </Typography>

          <Button
            disableElevation
            variant="contained"
            sizeVariant="large"
            sx={{
              color: 'white',
              background: 'black'
            }}
          >
            <Typography variant="pBodyBoldWt">
              Browse all projects
            </Typography>
          </Button>
        </Grid>
        <Box
          sx={{
            position: 'absolute'
            top: 0,
            left: 0
          }}>
          <Image
              src={projectExplorerBackshape}
              alt='Abstract plant'
          />
        </Box>
      </Box>
    </Box>
}
