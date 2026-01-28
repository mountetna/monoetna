'use client'

import * as React from 'react';
import Box from '@mui/system/Box';
import Grid from '@mui/system/Grid';
import Button from '@mui/material/Button';
import { useTheme }  from '@mui/material';
import Typography from '@mui/material/Typography';
import Image from '../../components/image/image';
import projectHeroImg from '/public/images/projects-hero-img.svg'

import FeatureProject from './feature-project';

export default function ProjectInfo() {
  const theme = useTheme();

  return <Box sx={{ padding: '90px 144px', height: 876 }}>
      <Typography
          variant='h2'
      >
          Projects
      </Typography>
      <Grid container sx={{ position: 'relative', marginTop: '25px' }}>
        <Grid size={6} item container direction='column' spacing={4}>
          <Grid item>
            <Typography variant="pLarge">
                The Data Library hosts a series of research projects containing multi-omic data and de-identified clinical data. Project data is shared in our research community to best leverage these rich datasets for biological discovery.
            </Typography>
          </Grid>

          <Grid item>
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
        </Grid>
        <Grid size={6} item container spacing={2}>
          <FeatureProject name='xhlt2' color={ theme.palette.orange.grade75 } />
          <FeatureProject name='ipi' color={ theme.palette.blue.grade75 } />
        </Grid>
        <Box
          sx={{
            position: 'absolute',
            top: -290,
            left: -95,
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
