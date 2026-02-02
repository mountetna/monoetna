import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import IconButton from '@mui/material/IconButton';
import Image from 'next/image';
import arrowUpRightLight from '/public/images/icons/arrow-up-right-light.svg';
import LinkoutButton from '../link/linkout-button.tsx';
import { useTheme } from '@mui/material';

import { useWindowDimensions } from '@/lib/utils/responsive';
import { ProjectExplorerContext } from './context';

import triangleDarkUp from '/public/images/icons/indicator-arrow-dark.svg'
import triangleLightUp from '/public/images/icons/indicator-arrow-light.svg'


export default function FeatureProject({
  name,
  color
}: {
  name: string;
  color: string;
}) {
    const { state: { projectData } } = React.useContext(ProjectExplorerContext);
    const project = projectData.find(p => p.name == name);
    return (
        <Box
            sx={{
                display: 'flex',
                flexDirection: 'column',
                justifyContent: 'space-between',
                p: '16px',
                borderRadius: '30px',
                color: 'unset',
                background: color
            }}
        >
            <Box
                sx={{
                  margin: '5px 0px'
                }}
            >
              <Typography sx={{ color: 'white' }} variant='h5'>
                  {project.fullName}
              </Typography>
            </Box>
            <Box
                sx={{
                    margin: '20px 20.5px',
                    borderRadius: '249px',
                    width: '249px',
                    height: '249px',
                    background: `url(${project.theme.coverImage.src})`,
                    backgroundSize: '363px',
                    backgroundPosition: '-57px -91px'
                }}
            >
            </Box>
            <Box
                sx={{
                  display: 'flex',
                  justifyContent: 'flex-end',
                  width: '100%'
                }}
            >
              <LinkoutButton size='large'
                tooltip='Open in Library'
                link={project.href}/>
            </Box>
        </Box>
    )
}
