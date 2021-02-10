import React, {useState, useEffect} from 'react';
import { isEditor } from 'etna-js/utils/janus';
import Icon from 'etna-js/components/icon';

const projectRoleKey = (p) => p.role + p.project_name_full.toUpperCase()

const UserProjects = ({projects, user}) => {
  console.log({projects});
  return <div id='user-projects'>
    <div className='title'>Your Projects</div>
      <div className='project header'>
        <div className='project_name'> Project Name </div>
        <div className='full_name'> Title </div>
        <div className='role'> Role </div>
        <div className='privileged' title='can see restricted data'> Privileged </div>
      </div>
    { projects.sort((p1,p2) => projectRoleKey(p1).localeCompare(projectRoleKey(p2))).map( project =>
      <div key={project.project_name} className='project'>
        <div className='project_name'>
        { isEditor(user, project.project_name) ?
          <a href={`/${project.project_name}`} >
            { project.project_name }
          </a> :
            <span>{ project.project_name }</span>
        }
        </div>
        <div className='full_name'>{ project.project_name_full }</div>
        <div className='role'>{ project.role }</div>
        <div className='privileged'>{ project.privileged && <Icon icon='check'/>}</div>
      </div>)
    }
  </div>
}

export default UserProjects;
