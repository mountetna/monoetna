import * as React from 'react';
import { getData } from '@/lib/clients/vesta-api/request';
import { Box, Typography } from '@mui/material';
import ProjectAccess from '@/components/project-access';

export default async function Access({params}:{
  params: {
    project_name: string;
  }
}) {
  const { projects } = await getData();
  const project = projects.find(p => p.name == params.project_name)
  

  if (!project) {
    return <Box sx={{px: 150, py: 150}}><Typography>Project {params.project_name} is not found.</Typography></Box>
  }

  return <ProjectAccess project={project} loginUrl={process.env.JANUS_URL} accessUrl={process.env.TIMUR_URL}/>
}
