'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Image from 'next/image';
import Typography from '@mui/material/Typography';
import Button from '@mui/material/Button';
import { projectDataTypes } from '@/lib/utils/filters';
import DataFieldChip from '@/components/data/data-field-chip';
import ThemeChip from '@/components/themes/theme-chip';
import ProjectPI from '@/components/project-explorer/project-pi'; 
import { Project, PrincipalInvestigator } from './project-explorer/models';
import libraryIcon from '/public/images/icons/library.svg'
import Link from './link/link';
import { useUser } from './user/context';

const ProjectViewer = ({project}:{ project: Project}) => {
  const dataTypes = projectDataTypes(project);

  const user = useUser();

  const hasAccess = user && (project.name in user.permissions);
  return <Box sx={{
    px: '144px',
    py: '50px'
  }}>
    <Box sx={{display: 'flex'}}> <ThemeChip theme={project.theme}/></Box>
    <Box sx={{ display: 'flex', alignItems: 'center', gap: '25px', py: '25px' }}>
      <Typography variant="h2">{project.fullName}</Typography>
      <DataFieldChip value={project.name.toUpperCase()}/>
    </Box>
    <Box sx={{ display: 'flex', mb: '45px', justifyContent: 'space-between' }}>
      <Box sx={{ width: '660px' }}>
        <Box sx={{display: 'flex', width: '100%' }}>
          <Box sx={{ display: 'flex', flex: '1 1 auto', flexDirection: 'column' }}>
            <Typography variant="p2XSBoldWt">PROJECT TYPE</Typography>
            <Typography variant="pMedium">{project.type}</Typography>
          </Box>
          <Box sx={{ display: 'flex', flex: '1 1 auto', flexDirection: 'column' }}>
            <Typography variant="p2XSBoldWt">START DATE</Typography>
            <Typography variant="pMedium">{project.startDate.getFullYear()}</Typography>
          </Box>
        </Box>
        <Box sx={{ py: '25px' }}><Typography variant="pMedium">{project.description}</Typography></Box>
        <Box sx={{ py: '15px' }}>
          <Button
            disableElevation
            variant="contained"
            size="large"
            sx={{
              color: 'white',
              bgcolor: hasAccess ? 'black' : 'orange.grade10',

              padding: '16px 32px'
            }}
            href={ project.href }
          >
          {
            hasAccess ?
            <>
              <Image
                style={{
                  filter: 'brightness(50)',
                  marginRight: '5px'
                }}
                width='24'
                height='24'
                src={libraryIcon}
                alt='open book'
              />
              <Typography variant="pBodyBoldWt">
                Open in Library
              </Typography>
            </>
            : <Link
                href={`/project/${project.name}/access`}
                sx={{
                  color: 'white'
                }}
            >
              <Typography variant="pBodyBoldWt">
                Get Access
              </Typography>
            </Link>
          }
          </Button>
        </Box>
      </Box>
      <Box sx={{ width: '660px' }}>
            <Image
              style={{
                marginRight: '5px',
                width: '660px',
                height: '488px',
                objectFit: 'cover',
                borderRadius: '24px'
              }}
              src={project.theme.coverImage}
              alt='open book'
            />
      </Box>
    </Box>

    <Box sx={{ position: 'relative', marginTop: '30px' }}>
      <Box sx={theme => ({
        background: theme.palette.orange.grade100,
        borderRadius: '8px',
        position: 'absolute',
        top: '-20px',
        mx: '60px',
        width: 'calc(100% - 48px - 75px)',
        clipPath: 'polygon(0 0, 100% 0, 100% 10px, 0 10px)',
        padding: '24px'
      })}>
      </Box>
      <Box sx={theme => ({
        background: theme.palette.orange.grade75,
        borderRadius: '8px',
        clipPath: 'polygon(0 0, 100% 0, 100% 10px, 0 10px)',
        position: 'absolute',
        width: 'calc(100% - 48px - 15px)',
        top: '-10px',
        mx: '30px',
        padding: '24px'
      })}>
      </Box>
      <Box sx={theme => ({
        background: theme.palette.orange.grade50,
        borderRadius: '16px',
        padding: '24px'
      })}>
        <Box sx={{
          display: 'flex'
        }}>
          <Box sx={{ width: '50%'}} >
            <Typography variant="h6BoldWt" color="white">Data Types</Typography>
            <Box sx={{
              display: 'flex',
              alignItems: 'center',
              flexWrap: 'wrap',
              gap: '10px',
              py: '15px'
            }}>
            {
              dataTypes.map( dt => <DataFieldChip key={dt} value={dt} color='utilityWhite.main'/> )
            }
            </Box>
          </Box>
          <Box sx={{ width: '50%'}} >
            <Typography variant="h6BoldWt" color="white">Principal Investigators</Typography>
            {
              project.principalInvestigators.map(
                (pi:PrincipalInvestigator) => <ProjectPI variant='filled' key={pi.name} data={pi}/>
              )
            }
          </Box>
        </Box>
        <Box sx={{ display: 'flex', flexDirection: 'column' }}> 
          <Typography variant="h6BoldWt" color="white">Project Status</Typography>
          <Typography variant="pMediumBoldWt" color="white">{project.status}</Typography>
        </Box>
      </Box>
    </Box>
  </Box>
}

export default ProjectViewer;
