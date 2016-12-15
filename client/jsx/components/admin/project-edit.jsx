import * as React from 'react';
import ProjectEntry from './project-entry';

export default class ProjectEdit extends React.Component{

  constructor(){

    super();
  }

  render(){

    var projects = this['props']['appState']['projects'];
    return (

      <div className='admin-edit-grouping'>

        <table id='project-edit-group'>

          <thead>

            <tr className='admin-edit-head-group'>

              <th id='project-name-column' className='admin-edit-title'>

                { 'project' }
                <div className='admin-column-head-arrow-group'>

                  <span className="glyphicon glyphicon-triangle-bottom"></span>
                </div>
              </th>
              <th id='project-group-column' className='admin-edit-title'>

                { 'group' }
                <div className='admin-column-head-arrow-group'>

                  <span className="glyphicon glyphicon-triangle-bottom"></span>
                </div>
              </th>
              <th id='project-control-column' className='admin-edit-title'>

                <button className='admin-add-btn'>

                  <span className='glyphicon glyphicon-plus white-glyphicon'></span>
                  { 'ADD PROJECT' }
                </button>
              </th>
            </tr>
          </thead>
          <tbody id='project-edit-body-group'>

            { projects.map((project, index)=>{

              return <ProjectEntry project={ project } key={ 'project-entry-'+ index } />
            }) }
          </tbody>
        </table>
      </div>
    );
  }
}