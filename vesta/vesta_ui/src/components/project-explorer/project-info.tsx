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

  return <Box sx={{ padding: '90px 144px', paddingBottom: '150px', background: `url(${projectHeroImg.src})`, backgroundPosition: '-100px -100px' }}>
      <Typography
          variant='h2'
      >
          Projects
      </Typography>
      <Box sx={{ display: 'flex', gap: '39px', position: 'relative', my: '25px', flexWrap: 'wrap' }}>
        <Box sx={{ display: 'flex', flexDirection: 'column', width: '600px', gap: '15px' }}>
          <Box>
            <Typography variant="pLarge">
                The Data Library hosts a series of research projects containing multi-omic data and de-identified clinical data. Project data is shared in our research community to best leverage these rich datasets for biological discovery.
            </Typography>
          </Box>
          <Box sx={{ py: '15px' }}>
            <Button
              href='#project-explorer'
              disableElevation
              variant="contained"
              size="large"
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
          </Box>
        </Box>
        <Box sx={{ width: '660px', display: 'flex', flexWrap: 'nowrap', gap: '16px' }}>
          <FeatureProject name='xhlt2' color={ theme.palette.orange.grade75 } />
          <FeatureProject name='ipi' color={ theme.palette.blue.grade75 } />
        </Box>
      </Box>
    </Box>
}
