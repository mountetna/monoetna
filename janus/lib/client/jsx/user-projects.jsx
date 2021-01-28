import React, {useState, useEffect} from 'react';
import { isEditor } from 'etna-js/utils/janus';

const projectRoleKey = (p) => p.role + p.project_name_full.toUpperCase()

const UserProjects = ({projects, user}) => {
  console.log({projects});
  return <div id='projects-group'>
    <div className='title'>Your Projects</div>
    { Object.values(projects).sort((p1,p2) => projectRoleKey(p1).localeCompare(projectRoleKey(p2))).map( project =>
      <div key={project.project_name} className='item'>
        { isEditor(user, project.project_name) ?
          <a href={`/${project.project_name}`} >
            { project.project_name_full }
          </a> :
            <span>{ project.project_name_full }</span>
        }
        - <i>{ project.role }{ project.privileged && ', privileged access'  }</i>
      </div>)
    }
  </div>
}

export default UserProjects;
