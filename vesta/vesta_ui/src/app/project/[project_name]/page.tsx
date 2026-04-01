import * as React from 'react';
import { getData } from '@/lib/clients/vesta-api/request';
import { Box, Typography } from '@mui/material';
import ProjectViewer from '@/components/project-viewer';

export default async function Project({params}:{
  params: {
    project_name: string;
  }
}) {
  const { projects } = await getData();
  const project = projects.find(p => p.name == params.project_name)

  if (!project) return <Box style={{ width: '100%', padding: '100px 150px' }}>
    <Typography variant="pMedium">No project named <b>{params.project_name}</b> is found!</Typography>
  </Box>;

  return <ProjectViewer project={project}/>
}
