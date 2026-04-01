'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import Link from './link/link';
import { useUser } from './user/context';
import { useRouter } from 'next/navigation';
import { Project } from './project-explorer/models';

const ProjectAccess = ({project, loginUrl, accessUrl}:{
  loginUrl: string;
  accessUrl: string;
  project: Project;
}) => {
  const user = useUser()

  const pi = project.principalInvestigators[0];
  const router = useRouter()

  const janusUrl = loginUrl + `/login?refer=${location.href}`;
  const timurUrl = accessUrl + '/' + project.name;

  if (user && (project.name in user.permissions || project.status == "Community" || project.status == "Resource")) {
    router.push(timurUrl, { scroll: true });
    return null;
  }

  return <Box
    sx={{
      px: '300px',
      py: '50px'
    }}
  >
    <Typography variant="h2">Get access to {project.name.toUpperCase()}</Typography>
    
    <Box
      sx={{
        p: '50px'
      }}>
        { !user && <ol>
          <li>
            <Typography variant="pMedium">If you already have a Library Card OR you have valid UCSF MyAccess credentials, login <Link href={ janusUrl }>here</Link>.</Typography>
          </li>
          <li><Typography variant="pMedium">
            If you do not have UCSF MyAccess credentials, contact the {project.name.toUpperCase()} principal investigator to request access to the project: <Link href={`mailto:${pi.email}?subject=Request access to ${project.name}`}>{pi.email}</Link>
            </Typography>
          </li>
        </ol>
        }
        { user &&
            <Typography variant="pMedium">
              You have a Library Card, but you don't have access to this project. Contact the {project.name.toUpperCase()} principal investigator to request access to the project: <Link href={`mailto:${pi.email}?subject=Request access to ${project.name}`}>{pi.email}</Link>
            </Typography>
        }
    </Box>
  </Box>
}

export default ProjectAccess;
