import * as React from 'react';
import { getData } from '@/lib/clients/vesta-api/request';
import { Box, Typography } from '@mui/material';
import ProjectAccess from '@/components/project-access';

type Params = Promise<{ project_name: string }>;

export default async function Access({params}:{
  params: Params 
}) {
  const { projects } = await getData();
  const { project_name } = await params;
  const project = projects.find(p => p.name == project_name)
  

  if (!project) {
    return <Box sx={{px: 150, py: 150}}><Typography>Project {project_name} is not found.</Typography></Box>
  }

  return <ProjectAccess project={project} loginUrl={process.env.JANUS_URL} accessUrl={process.env.TIMUR_URL}/>
}
