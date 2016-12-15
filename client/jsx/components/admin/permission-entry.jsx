import * as React from 'react';
import GenericAdminEntry from './generic-admin-entry';
import ProjectSearchDropdown from './project-search-dropdown';
import PermissionDropdown from './permission-dropdown';

export default class PermissionEntry extends GenericAdminEntry{

  componentDidMount(){

    var permission = this['props']['permission'];
    if(String(permission['id']).indexOf('permission') !== -1){

      this.setState({

        'editShow': true,
        'editActive': true,
        'permission': this['props']['permission']
      });
    }
  }

  checkForPrimary(){

    var perm = this['props']['permission'];
    if(perm['role'] == 'administrator'){

      if(perm['projectName'] == "administration"){

        alert('You cannot edit the primary permission.');
        return false;
      }
    }

    return true;
  }

  activateEntryEdit(){

    if(!this.checkForPrimary()) return;
    this.setState({ 'editActive': true });
  }

  resetEntry(){

    if(!this.checkForPrimary()) return;
  }

  deactivateEntryEdit(){

    this.setState({ 'editActive': false });

    var permission = this['props']['permission'];    
    if(permission['id'].indexOf('permission') !== -1){

      //remove the new and unsaved permission
      this['props']['callbacks'].removeUnsavedPermission(permission['id']);
    }
  }

  saveEntry(){

    if(!this.checkForPrimary()) return;
  }

  deleteEntry(){

    if(!this.checkForPrimary()) return;

    var permission = this['props']['permission'];
    if(permission['id'].indexOf('permission') !== -1){

      //remove the new and unsaved permission
      this['props']['callbacks'].removeUnsavedPermission(permission['id']);
    }
  }

  updateFromInput(field, value){

    //console.log(field, value);
    //var permission = this['props']['permission'];
    //this['props']['callbacks'].update
    //var permission = this['props']['permission'];

    //var permission = this['state']['permission'];
    //if(permission != undefined){
      
    //  permission[field] = value;
    //  this.setState({ 'permission': permission });
    //}
  }

  render(){

    var permission = this['props']['permission'];

    var callbacks = {

      'updateFromInput': this['updateFromInput'].bind(this)
    };

    var adminEntryProps = {

      'className': 'admin-edit-entry-group',
      'onMouseEnter': this['showControlGroup'].bind(this),
      'onMouseLeave': this['hideControlGroup'].bind(this)
    };

    var userEmailEntryProps = {

      'className': 'admin-entry-input',
      'defaultValue': this['props']['permission']['userEmail']
    };

    if(!this['state']['editActive']){

      userEmailEntryProps['className'] = 'admin-entry-input-inactive';
      userEmailEntryProps['disabled'] = true;
    }

    var projectDropdownProps = {

      'permission': this['props']['permission'],
      'projects': this['props']['projects'],
      'editActive': this['state']['editActive'],
      'callbacks': callbacks
    };

    var permDropdownProps = {

      'permission': this['props']['permission'],
      'editActive': this['state']['editActive'],
      'callbacks': callbacks
    };

    return (

      <tr { ...adminEntryProps }>

        <td>

          <input { ...userEmailEntryProps } />
        </td>
        <td>

          <ProjectSearchDropdown { ...projectDropdownProps } />
        </td>
        <td>

          <PermissionDropdown { ...permDropdownProps } />
        </td>
        <td>

          { this.renderEditControlGroup() }
        </td>
      </tr>
    );
  }
}