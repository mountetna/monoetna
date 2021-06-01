import React, {useState, useEffect, useCallback} from 'react';
import Icon from 'etna-js/components/icon';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {selectUser} from 'etna-js/selectors/user-selector';
import { json_post, json_get } from 'etna-js/utils/fetch';
import {isSuperuser} from 'etna-js/utils/janus';

import TextField from '@material-ui/core/TextField';
import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableContainer from '@material-ui/core/TableContainer';
import TableHead from '@material-ui/core/TableHead';
import TableRow from '@material-ui/core/TableRow';
import Paper from '@material-ui/core/Paper';

const Projects = ({projects}) => (
  <div id='admin-projects'>
    <div className='title'>Projects</div>
    <TableContainer component={Paper}>
      <Table aria-label='all projects'>
        <TableHead>
          <TableRow>
            <TableCell>Project Name</TableCell>
            <TableCell>Title</TableCell>
          </TableRow>
        </TableHead>
        <TableBody>
          { projects.sort((p1,p2) => p1.project_name_full.localeCompare(p2.project_name_full)).map( project =>
            <TableRow key={project.project_name}>
              <TableCell>
                <a href={`/${project.project_name}`} >
                  { project.project_name }
                </a>
              </TableCell>
              <TableCell>{ project.project_name_full }</TableCell>
            </TableRow>)
          }
        </TableBody>
      </Table>
    </TableContainer>
  </div>
);

const postAddProject = (project) => json_post('/add_project', project);

const NewProject = ({retrieveAllProjects}) => {
  let [ newproject, setNewProject ] = useState({})
  let [ error, setError ] = useState(null);

  return <div id='new-project'>
    <div className='title'>New Project</div>
    { error && <div className='error'>Error: {error}</div> }
    <div className='item'>
      <div className='cell'>
        <TextField
          placeholder='Project Full Name'
          value={ newproject.project_name_full || '' }
          onChange={ (e) => setNewProject({ ...newproject, project_name_full: e.target.value }) }/>
      </div>
      <div className='cell'>
        <TextField placeholder='project_short_name'
          value={ newproject.project_name || '' }
          onChange={ (e) => setNewProject({ ...newproject, project_name: e.target.value }) }/>
      </div>
      <div className='cell submit'>
        <Icon className='approve' icon='magic' onClick={
          () => postAddProject( newproject ).then(
            ()=> {
              retrieveAllProjects();
              setNewProject({});
              setError(null);
            }
          ).catch(
            e => e.then( ({error}) => setError(error) )
          )
        }/>
      </div>
    </div>
  </div>
}

const JanusAdmin = () => {
  let user = useReduxState( state => selectUser(state) );

  let [ projects, setProjects ]  = useState([]);

  let retrieveAllProjects = useCallback(
    () => {
      json_get('/allprojects')
        .then(({projects}) => setProjects(projects))
    }, [ setProjects ]
  );

  useEffect(retrieveAllProjects, []);
  return <div id='janus-admin'>
    <Projects user={user} projects={projects}/>
    { isSuperuser(user) && <NewProject retrieveAllProjects={retrieveAllProjects}/> }
  </div>;
}

export default JanusAdmin;
