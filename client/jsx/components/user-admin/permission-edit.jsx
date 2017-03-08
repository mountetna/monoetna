import * as React from 'react';
import PermissionEntry from './permission-entry';

export default class PermissionEdit extends React.Component{

  constructor(){

    super();
  }

  selectFile(){

    /*
     * We are using a button to surragate the file input so we may have 
     * a custom browse button.
     */
    document.getElementById('permission-file-selector').click();
  }

  fileSelected(event){

    if(event === undefined) return;
    var fileSelector = event.target;
    var file = fileSelector.files[0];
    this['props'].uploadPermissions(file);
  }

  render(){

    var permissions = this['props']['adminInfo']['permissions'];
    var projects = this['props']['adminInfo']['projects'];
    var users = this['props']['adminInfo']['users'];

    var permFileSelector = {

      'id': 'permission-file-selector',
      'type': 'file',
      'name': 'upload-file',
      'onChange': this['fileSelected'].bind(this)
    };

    var uploadPermBtnProps = {

      'className': 'admin-add-btn',
      'onClick': this['selectFile'].bind(this)
    };

    var dwnldPermBtnProps = {

      'className': 'admin-add-btn',
      'onClick': this['props']['downloadPermissions'].bind(this)
    };

    var addPermBtnProps = {

      'className': 'admin-add-btn',
      'onClick': this['props']['addPermission'].bind(this)
    };

    var callbacks = {

      'removeUnsavedPermission': this['props']['removeUnsavedPermission'].bind(this),
      'savePermission': this['props']['savePermission'].bind(this)
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
                  { ' ADD PERM' }
                </button>
                <button { ...dwnldPermBtnProps }>

                  <span className='glyphicon glyphicon-save-file white-glyphicon'></span>
                  { ' DOWN' }
                </button>
                <button { ...uploadPermBtnProps }>

                  <input { ...permFileSelector } />
                  <span className='glyphicon glyphicon-open-file white-glyphicon'></span>
                  { ' UP' }
                </button>
              </th>
            </tr>
          </thead>
          <tbody id='permission-edit-body-group'>

            { permissions.map((permission, index)=>{

              var permEntryProps = {

                'permission': permission,
                'projects': projects,
                'users': users,
                'key': permission['reactKey'],
                'reactKey': permission['reactKey'],
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