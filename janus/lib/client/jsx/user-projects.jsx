import React, {useState, useEffect} from 'react';
import { isEditor } from 'etna-js/utils/janus';
import Icon from 'etna-js/components/icon';

import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableContainer from '@material-ui/core/TableContainer';
import TableHead from '@material-ui/core/TableHead';
import TableRow from '@material-ui/core/TableRow';
import Paper from '@material-ui/core/Paper';

const projectRoleKey = (p) => p.role + p.project_name_full.toUpperCase()

const UserProjects = ({projects, user}) => {
  return <div id='user-projects'>
    <div className='title'>Your Projects</div>
    <TableContainer component={Paper}>
      <Table aria-label='project users'>
        <TableHead>
          <TableRow>
            <TableCell>Project Name</TableCell>
            <TableCell>Title</TableCell>
            <TableCell>Role</TableCell>
            <TableCell padding='checkbox'>Privileged</TableCell>
          </TableRow>
        </TableHead>
        <TableBody>
          { projects.sort((p1,p2) => projectRoleKey(p1).localeCompare(projectRoleKey(p2))).map( project =>
            <TableRow key={project.project_name}>
              <TableCell>
              { isEditor(user, project.project_name) ?
                <a href={`/${project.project_name}`} >
                  { project.project_name }
                </a> : project.project_name
              }
              </TableCell>
              <TableCell>{ project.project_name_full }</TableCell>
              <TableCell>{ project.role }</TableCell>
              <TableCell className='privileged'>{ project.privileged && <Icon icon='check'/>}</TableCell>
            </TableRow>)
          }
        </TableBody>
      </Table>
    </TableContainer>
  </div>
}

export default UserProjects;
