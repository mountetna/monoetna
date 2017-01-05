import * as React from 'react';
import GenericAdminEntry from './generic-admin-entry';
import ProjectSearchDropdown from './project-search-dropdown';
import PermissionDropdown from './permission-dropdown';

export default class PermissionEntry extends GenericAdminEntry{

  componentDidMount(){

    var permission = this['props']['permission'];
    if(permission['projectId'] == null){

      this.setState({

        'editShow': true,
        'editActive': true
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

  deactivateEntryEdit(){

    var permission = this['props']['permission'];
    if(permission['id'] == null) return;

    this.setState({ 'editActive': false });
  }

  resetEntry(){

    if(!this.checkForPrimary()) return;
    //this.forceUpdate();
  }

  saveEntry(){

    if(!this.checkForPrimary()) return;

    var reactKey = this['props']['permission']['reactKey'];
    var email = this['userEmailInput']['value'];
    var prjNm = this['projectDropdownComponent']['state']['inputValue'];
    var role = this['permissionDropdownComponent']['state']['inputValue'];

    /*
     * These are only simple validations to keep the UI tidy. There are more 
     * stringent validations higher up in the UI and definately on the server.
     */
    if(!this.validateUser(email)){

      alert('Please select a valid user.');
      return;
    }

    if(!this.validateProject(prjNm)){

      alert('Please select a valid project.');
      return;
    }

    if(!this.validateRole(role)){

      alert('Please select a permission');
      return;
    }

    var permission = this['props']['permission'];
    permission['userEmail'] = email;
    permission['projectName'] = prjNm;
    permission['role'] = role;

    this['props']['callbacks'].savePermission(permission);
  }

  deleteEntry(){

    if(!this.checkForPrimary()) return;
    var reactKey = this['props']['permission']['reactKey'];
    this['props']['callbacks'].removeUnsavedPermission(reactKey);
  }

  validateUser(email){

    if(email == '') return false;
    if(!VALIDATE_EMAIL(email)) return false;

    return true;
  }

  validateProject(projectName){

    if(projectName == '') return false;

    var exists = false;
    for(var a = 0; a < this['props']['projects']['length']; ++a){

      var prjNm = this['props']['projects'][a]['projectName'].toLowerCase();
      if(projectName == prjNm){

        exists = true;
      }
    }

    return exists;
  }

  validateRole(role){

    if(role == '') return false;

    return true;
  }

  render(){

    var permission = this['props']['permission'];
    var adminEntryProps = {

      'className': 'admin-edit-entry-group',
      'onMouseEnter': this['showControlGroup'].bind(this),
      'onMouseLeave': this['hideControlGroup'].bind(this)
    };

    var userEmailEntryProps = {

      'className': 'admin-entry-input',
      'defaultValue': this['props']['permission']['userEmail'],
      'ref': (input)=>{ this['userEmailInput'] = input }
    };

    if(!this['state']['editActive']){

      userEmailEntryProps['className'] = 'admin-entry-input-inactive';
      userEmailEntryProps['disabled'] = true;
    }

    var projectDropdownProps = {

      'permission': this['props']['permission'],
      'projects': this['props']['projects'],
      'editActive': this['state']['editActive'],
      'ref': (component)=>{ this['projectDropdownComponent'] = component }
    };

    var permDropdownProps = {

      'permission': this['props']['permission'],
      'editActive': this['state']['editActive'],
      'ref': (component)=>{ this['permissionDropdownComponent'] = component }
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