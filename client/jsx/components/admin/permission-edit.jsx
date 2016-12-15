import * as React from 'react';
import PermissionEntry from './permission-entry';

export default class PermissionEdit extends React.Component{

  constructor(){

    super();
  }

  addNewPermission(){

    this['props'].addPermission();
  }

  removeUnsavedPermission(permissionId){

    this['props'].removeUnsavedPermission(permissionId);
  }

  render(){

    var permissions = this['props']['appState']['permissions'];
    var projects = this['props']['appState']['projects'];

    var addPermBtnProps = {

      'className': 'admin-add-btn',
      'onClick': this['addNewPermission'].bind(this)
    };

    var callbacks = {

      'removeUnsavedPermission': this['removeUnsavedPermission'].bind(this)
    };
    
    return (

      <div className='admin-edit-grouping'>

        <table id='permission-edit-group'>

          <thead>

            <tr className='admin-edit-head-group'>

              <th id='permission-user-column' className='admin-edit-title'>

                { 'user email ' }
                <div className='admin-column-head-arrow-group'>

                  <span className="glyphicon glyphicon-triangle-bottom"></span>
                </div>
              </th>
              <th id='permission-project-column' className='admin-edit-title'>

                { 'project ' }
                <div className='admin-column-head-arrow-group'>

                  <span className="glyphicon glyphicon-triangle-bottom"></span>
                </div>
              </th>
              <th id='permission-column' className='admin-edit-title'>

                { 'permission ' }
                <div className='admin-column-head-arrow-group'>

                  <span className="glyphicon glyphicon-triangle-bottom"></span>
                </div>
              </th>
              <th id='permission-control-column' className='admin-edit-title'>

                <button { ...addPermBtnProps }>

                  <span className='glyphicon glyphicon-plus white-glyphicon'></span>
                  { 'ADD PERMISSION' }
                </button>
              </th>
            </tr>
          </thead>
          <tbody id='permission-edit-body-group'>

            { permissions.map((permission, index)=>{

              var permEntryProps = {

                'permission': permission,
                'projects': projects,
                'key': permission['reactKey'],
                'callbacks': callbacks
              };

              return <PermissionEntry { ...permEntryProps } />
            }) }
          </tbody>
        </table>
      </div>
    );
  }
}