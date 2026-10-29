import * as React from 'react';
import { getData } from '@/lib/clients/vesta-api/request';
import { Box, Typography } from '@mui/material';
import ProjectViewer from '@/components/project-viewer';

type Params = Promise<{ project_name: string }>;

export default async function Project({params}:{
  params: Params
}) {
  const { projects } = await getData();
  const { project_name } = await params;
  const project = projects.find(p => p.name == project_name)

  if (!project) return <Box style={{ width: '100%', padding: '100px 150px' }}>
    <Typography variant="pMedium">No project named <b>{project_name}</b> is found!</Typography>
  </Box>;

  return <ProjectViewer project={project}/>
}
