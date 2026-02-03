'use client'

import * as React from 'react';
import { ProjectExplorerContext } from './context';
import ProjectHeadingInfo from './project-heading-info';
import {Project, ProjectHeadingInfoSet } from '@/components/project-explorer/models';

import Box from '@mui/system/Box'
import Table from '@mui/material/Table';
import TableBody from '@mui/material/TableBody';
import TableCell from '@mui/material/TableCell';
import TableContainer from '@mui/material/TableContainer';
import TableHead from '@mui/material/TableHead';
import TableRow from '@mui/material/TableRow';
import Paper from '@mui/material/Paper';
import Pagination, { PaginationClasses } from '@/components/searchable-list/controls/pagination'
import ThemeChip from '@/components/themes/theme-chip';
import DataFieldChip from '@/components/data/data-field-chip';
import Typography from '@mui/material/Typography';
import { useTheme } from '@mui/material';
import { projectDataTypes } from '@/lib/utils/filters';
import LinkoutButton from '../link/linkout-button';

const SortTriangle = ({selected, onClick}:{
  selected: boolean,
  onClick: React.MouseEventHandler<HTMLDivElement>
}) => {
  const theme = useTheme();
  return <Box sx={{ display: 'inline', px: '4px', cursor: 'pointer' }} onClick={ onClick }>
    <svg width="17" height="9" viewBox="0 0 17 9" fill="none" xmlns="http://www.w3.org/2000/svg">
      <path d="M8.5 9L17 0H0L8.5 9Z" fill={ selected ? theme.palette.ground.grade50 : theme.palette.ground.grade75 }/>
    </svg>
  </Box>
}

const DataTypeColumn = ({project}:{project:Project}) => {
  const extra = 3;
  const dataTypes = projectDataTypes(project);

  const allTypes = dataTypes.sort();

  const shownTypes = allTypes.slice(0, extra);

  return <Box
  sx={{ minWidth: '600px',
  columnGap: '10px',
  display: 'flex'}}>
  {
    shownTypes.map(
      dataType => <DataFieldChip key={dataType} value={dataType} />
    )
  }
  {
    allTypes.length > extra ? <DataFieldChip value={`+${ allTypes.length - extra }`}/> : null
  }
  </Box>
}

const columnComponent = (columnName:string, project:Project) => {
  switch (columnName) {
    case "Project ID":
      return <DataFieldChip value={project.name.toUpperCase()}/>;
    case "Name":
      return <Typography sx={{
        whiteSpace: 'nowrap',
        maxWidth: '250px',
        overflow: 'hidden',
        textOverflow: 'ellipsis'

      }} variant="pBodyMediumWt">{project.fullName}</Typography>;
    case "Data types":
      return <DataTypeColumn project={project}/>;
    case "Investigators":
      return <ProjectHeadingInfo
        projectData={project}
        projectOpen={false}
        infoSet={ ProjectHeadingInfoSet.pis }
        variant='small'
      />
    case "Subjects":
      return <Typography variant="pBodyMediumWt">{project.subjectCount}</Typography>
    case "Year":
      return <Typography variant="pBodyMediumWt">{project.startDate.getFullYear()}</Typography>
    case "Type":
      return <Typography variant="pBodyMediumWt">{project.type}</Typography>
    case "Phase":
      return <Typography variant="pBodyMediumWt">{project.status}</Typography>
    case "Theme":
      return <ThemeChip theme={project.theme}/>
    default:
      return columnName;
  }
}


const ProjectTable = ({currentPage, setCurrentPage}:{
  currentPage: number,
  setCurrentPage: (newPage:number) => void
}) => {
  const { state: { projectData, visibleColumns }, filteredProjectData } = React.useContext(ProjectExplorerContext);
  const [sortColumn, setSortColumn] = React.useState('Name');
  const pageSize = 8;
  const listSize = filteredProjectData.length;
  const theme = useTheme();

  const handleSetCurrentPage = React.useCallback(
    (newPage:number) => newPage > -1 && newPage < listSize && setCurrentPage(newPage),
    [ currentPage, listSize, setCurrentPage ]
  );

  const projectsPage = filteredProjectData.slice(
    Math.min(currentPage * pageSize + 1, listSize) - 1,
    Math.min(currentPage * pageSize + pageSize, listSize + 1) - 1
  );

  const paginationEl = (
      <Pagination
          currentPage={currentPage}
          pageSize={pageSize}
          listSize={listSize}
          onClickPrev={() => handleSetCurrentPage(currentPage - 1)}
          onClickNext={() => handleSetCurrentPage(currentPage + 1)}
          listItemLabel='projects'
          sx={{
              [theme.breakpoints.up('tablet')]: {
                  height: '100%',
                  [`& .${PaginationClasses.pageInfo}`]: {
                      height: '100%',
                  },
              },
          }}
      />
  )

  return (
  <>
    <TableContainer sx={{ height: '600px' }} >
      <Table aria-label="project table">
        <TableHead>
          <TableRow sx={ theme => ({ borderBottom: `1px solid ${theme.palette.ground.grade50}`}) }>
            {
              visibleColumns.map(
                (columnName:string) =>
                <TableCell 
                  sx={{
                    bgcolor: 'utilityHighlight.main',
                    borderBottom: 'none',
                    minWidth: '130px'
                  }}
                  key={columnName}>
                  <Typography variant="pBodyMediumWt">{columnName}</Typography>
                  <SortTriangle
                    selected={ columnName == sortColumn }
                    onClick={ () => setSortColumn( columnName ) }/>
                </TableCell>
              )
            }
            <TableCell sx={{
              bgcolor: 'utilityHighlight.main',
              borderBottom: 'none',
            }}></TableCell>
          </TableRow>
        </TableHead>
        <TableBody>
        {
          projectsPage.map( (project:Project) =>
          <TableRow key={project.name} sx={{
                    height: '73px'
          }}>
            {
              visibleColumns.map(
                (columnName:string) => <TableCell
                  sx={{
                    py: '10px',
                    bgcolor: 'utilityHighlight.main',
                    borderBottom: 'none'
                  }}
                  key={columnName}>
                  <Box sx={{ display: 'flex' }}>
                  {columnComponent( columnName, project )}
                  </Box>
                </TableCell>
              )
            }
            <TableCell sx={{
              bgcolor: 'utilityHighlight.main',
              borderBottom: 'none'
            }}>
              <LinkoutButton size='small'
                tooltip='Open in Library'
                link={ project.href }
                />
            </TableCell>
          </TableRow>
          )
        }
        </TableBody>
      </Table>
    </TableContainer>
    <Box sx={{ width: '50%'}}>
      {paginationEl}
    </Box>
  </>
  );
}

export default ProjectTable;
