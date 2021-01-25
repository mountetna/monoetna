import React, {useState, useEffect} from 'react';

const projectRoleKey = (p) => p.role + p.project.project_name_full.toUpperCase()

const UserProjects = ({projects}) => {
  console.log({projects});
  <div id='projects-group'>
    <div class='title'>Your Projects</div>
    { Object.values(projects).sort((p1,p2) => projectRoleKey(p1).localeCompare(projectRoleKey(p2))).map( project =>
      <div key={project.project_name} class='item'>
        { project.editor ?
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
