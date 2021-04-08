import React, {useState, useEffect, useCallback} from 'react';
import Icon from 'etna-js/components/icon';
import {useReduxState} from 'etna-js/hooks/useReduxState';
import {selectUser} from 'etna-js/selectors/user-selector';
import { json_post, json_get } from 'etna-js/utils/fetch';
import {isSuperuser} from 'etna-js/utils/janus';

const Projects = ({projects}) => (
  <div id='admin-projects'>
    <div className='title'>Projects</div>
      <div className='project header'>
        <div className='project_name'> project_name </div>
        <div className='full_name'> title </div>
      </div>
    { projects.sort((p1,p2) => p1.project_name_full.localeCompare(p2.project_name_full)).map( project =>
      <div key={project.project_name} className='project'>
        <div className='project_name'>
          <a href={`/${project.project_name}`} >
            { project.project_name }
          </a>
        </div>
        <div className='full_name'>{ project.project_name_full }</div>
      </div>)
    }
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
        <input type='text' placeholder='Project Full Name' name='project_name_full'
          value={ newproject.project_name_full || '' }
          onChange={ (e) => setNewProject({ ...newproject, project_name_full: e.target.value }) }/>
      </div>
      <div className='cell'>
        <input type='text' placeholder='project_short_name' name='project_name'
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
