import * as React from 'react';
import { getData } from '@/lib/clients/vesta-api/request';
import { Box, Typography } from '@mui/material';
import ProjectAccess from '@/components/project-access';

export default async function Access({params}) {
  const { projects } = await getData();
  const project = projects.find(p => p.name == params.project_name)

  return <ProjectAccess project={project} loginUrl={process.env.JANUS_URL} accessUrl={process.env.TIMUR_URL}/>
}
