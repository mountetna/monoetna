import * as React from 'react';
import GenericAdminEntry from './generic-admin-entry';
import ProjectSearchDropdown from './project-search-dropdown';
import RoleDropdown from './role-dropdown';

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

  deleteEntry(){

    if(!this.checkForPrimary()) return;
    var reactKey = this['props']['permission']['reactKey'];
    this['props']['callbacks'].removeUnsavedPermission(reactKey);
  }

  saveEntry(){

    if(!this.checkForPrimary()) return;

    /*
     * These are only simple validations to keep the UI tidy. There are more 
     * stringent validations higher up in the UI and definately on the server.
     */
    var email = this['userEmailInput']['value'];
    var userId = this.validateUser(email);
    if(userId <= 0){

      alert('Please select a valid user.');
      return;
    }
    this['props']['permission']['userId'] = userId;
    this['props']['permission']['userEmail'] = email;

    var projectName = this['props']['permission']['projectName'];
    var projectId = this.validateProject(projectName);
    if(projectId <= 0){

      alert('Please select a valid project.');
      return;
    }

    var role = this['props']['permission']['role'];
    if(!this.validateRole(role)){

      alert('Please select a permission');
      return;
    }

    this['props']['callbacks'].savePermission(this['props']['permission']);
  }

  projectSelected(project){

    var projectId = this.validateProject(project);
    if(projectId <= 0){

      alert('Please select a valid project.');
      return;
    }

    this['props']['permission']['projectName'] = project;
    this['props']['permission']['projectId'] = projectId;
  }

  roleSelected(role){

    if(!this.validateRole(role)){

      alert('Please select a permission');
      return;
    }

    this['props']['permission']['role'] = role;
  }

  validateUser(email){

    if(email == '') return false;
    if(!VALIDATE_EMAIL(email)) return false;

    var userId = 0;
    for(var a = 0; a < this['props']['users']['length']; ++a){

      var userEmail = this['props']['users'][a]['email'].toLowerCase();
      if(email == userEmail){

        userId = this['props']['users'][a]['id'];
      }
    }

    return userId;
  }

  validateProject(projectName){

    if(projectName == '') return false;

    var projectId = 0;
    for(var a = 0; a < this['props']['projects']['length']; ++a){

      var prjNm = this['props']['projects'][a]['projectName'].toLowerCase();
      if(projectName == prjNm){

        projectId = this['props']['projects'][a]['id'];
      }
    }

    return projectId;
  }

  validateRole(role){

    if(role == '') return false;

    return true;
  }

  render(){

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
      'callbacks': {

        'projectSelected': this['projectSelected'].bind(this)
      }
    };

    var roleDropdownProps = {

      'permission': this['props']['permission'],
      'editActive': this['state']['editActive'],
      'callbacks': {

        'roleSelected': this['roleSelected'].bind(this)
      }
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

          <RoleDropdown { ...roleDropdownProps } />
        </td>
        <td>

          { this.renderEditControlGroup() }
        </td>
      </tr>
    );
  }
}