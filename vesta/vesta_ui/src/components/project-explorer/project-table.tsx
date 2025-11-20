'use client'

import * as React from 'react';
import { ProjectExplorerContext } from './context';

import Box from '@mui/system/Box'
import Table from '@mui/material/Table';
import TableBody from '@mui/material/TableBody';
import TableCell from '@mui/material/TableCell';
import TableContainer from '@mui/material/TableContainer';
import TableHead from '@mui/material/TableHead';
import TableRow from '@mui/material/TableRow';
import Paper from '@mui/material/Paper';
import Pagination, { PaginationClasses } from '@/components/searchable-list/controls/pagination'
import Typography from '@mui/material/Typography';
import { useTheme } from '@mui/material';

const columnComponent = (columnName, project) => {
  switch (columnName) {
    case "Name":
      return <Typography>{project.name}</Typography>
    default:
      return columnName;
  }
}


const ProjectTable = () => {
  const { state: { projectData, visibleColumns } } = React.useContext(ProjectExplorerContext);
  const [currentPage, setCurrentPage] = React.useState(0);
  const pageSize = 12;
  const listSize = projectData.length;
  const theme = useTheme();

  const handleSetCurrentPage = React.useCallback(
    newPage => newPage > -1 && newPage < listSize && setCurrentPage(newPage),
    [ currentPage, listSize, setCurrentPage ]
  );

  const projectsPage = projectData.slice(
    Math.min(currentPage * pageSize + 1, listSize) - 1,
    Math.min(currentPage * pageSize + pageSize, listSize) - 1
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
    <TableContainer component={Paper}>
      <Table aria-label="project table">
        <TableHead>
          <TableRow>
            {
              visibleColumns.map(
                columnName => <TableCell key={columnName}>{columnName}</TableCell>
              )
            }
          </TableRow>
        </TableHead>
        <TableBody>
        {
          projectsPage.map( project =>
          <TableRow>
            {
              visibleColumns.map(
                columnName => <TableCell key={columnName}>{columnComponent( columnName, project )}</TableCell>
              )
            }
          </TableRow>
          )
        }
        </TableBody>
      </Table>
      <Box>
        {paginationEl}
      </Box>
    </TableContainer>
  );
}

export default ProjectTable;
