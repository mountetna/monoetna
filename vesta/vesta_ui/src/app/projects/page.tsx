import * as React from 'react';
import Box from '@mui/system/Box';
import ProjectExplorer from '@/components/project-explorer/project-explorer';
import {ProjectExplorerContextProvider} from '@/components/project-explorer/context';
import ProjectInfo from '@/components/project-explorer/project-info';
import { getData } from '@/lib/clients/vesta-api/request';
import _ from 'lodash';

export default async function Projects() {
  const data = await getData();
  return (
    <Box>
      <ProjectExplorerContextProvider projectData={_.sortBy(data.projects, (p) => p.fullName)}>
        <ProjectInfo/>
        <ProjectExplorer/>
      </ProjectExplorerContextProvider>
    </Box>
  );
}
