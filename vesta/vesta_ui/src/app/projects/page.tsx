import * as React from 'react';
import Box from '@mui/system/Box';
import _ from 'lodash';
import { getData } from '@/lib/clients/vesta-api/request';
import ProjectExplorer from '@/components/project-explorer/project-explorer';
import ProjectInfo from '@/components/project-explorer/project-info';

export default async function Projects() {
  const data = await getData()

  return (
    <Box>
      <ProjectInfo/>

      <ProjectExplorer
        projectData={_.sortBy(data.projects, (p) => p.fullName)}
      />
    </Box>
  );
}
