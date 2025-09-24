'use client'

import * as React from 'react';
import Box from '@mui/system/Box';
import Grid from '@mui/system/Grid';
import Button from '@mui/material/Button';
import Typography from '@mui/material/Typography';
import Image from '../../components/image/image';
import projectHeroImg from '/public/images/projects-hero-img.svg'

export default function ProjectInfo() {
  return <Box sx={{ padding: '90px 144px', height: 876 }}>
      <Typography
          variant='h2'
      >
          Projects
      </Typography>
      <Grid container sx={{ position: 'relative', marginTop: '25px' }}>
        <Grid size={6} item container spacing={8}>
          <Typography variant="pLarge">
              The Data Library hosts a series of research projects containing multi-omic data and de-identified clinical data. Project data is shared in our research community to best leverage these rich datasets for biological discovery.
          </Typography>

          <Button
            disableElevation
            variant="contained"
            sizeVariant="large"
            sx={{
              color: 'white',
              background: 'black',
              padding: '16px 32px'
            }}
          >
            <Typography variant="pBodyBoldWt">
              Browse all projects
            </Typography>
          </Button>
        </Grid>
        <Box
          sx={{
            position: 'absolute',
            top: -200,
            left: -100,
            zIndex: -1
          }}>
          <Image
              src={projectHeroImg}
              alt='Abstract plant'
          />
        </Box>
      </Grid>
    </Box>
}
