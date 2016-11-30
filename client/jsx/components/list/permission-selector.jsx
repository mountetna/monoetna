import * as React from 'react';

export default class PermissionSelector extends React.Component{

  constructor(props){

    super(props);

    this['state'] = {

      'projectTrayActive': false,
      'inputValue': '',
      'projectRole': '',
    };
  }

  toggleTray(event = undefined){

    var trayActive = (this['state']['projectTrayActive']) ? false : true;
    this.setState({ 'projectTrayActive': trayActive });
  }

  disableTray(event = undefined){

    this.setState({ 'projectTrayActive': false });
  }

  trayEntrySelected(event = undefined){

    var projectName = event['target']['dataset']['name'];
    var projectId = event['target']['dataset']['id'];
    var role = event['target']['dataset']['role'];

    this.setState({ 

      'projectRole': role, 
      'inputValue': projectName,
      'projectTrayActive': false
    });

    // Bubble up the data to the parent.
    this['props']['callbacks'].projectSelected(projectName, role);
  }

  updateFileName(event){

    var value = event['target']['value'];
    var role = this.getRole(value);

    // If the search length is longer than 3 show the dropdown.
    if(value['length'] >= 3 && !this['state']['projectTrayActive']){

      this.setState({ 

        'inputValue': value,
        'projectTrayActive': true,
        'projectRole': role
      });
    }
    else{

      this.setState({ 

        'inputValue': value,
        'projectRole': role
      });
    }

    /*
     * If the role is not empty then we have a valid project and can update the
     * file data.
     */
    if(role != ''){

      // Bubble up the data to the parent.
      this['props']['callbacks'].projectSelected(projectName, role);
    }
  }

  getRole(projectName){

    var permissions = this['props']['permissions'];
    var role = '';
    for(var index in permissions){

      if(permissions[index]['projectName'] == projectName){

        role = permissions[index]['role'];
      }
    }

    return role;
  }

  filterPermissions(){

    // The list of entries to be 'searched' over.
    var permissions = this['props']['permissions'];

    // The value of the input to 'search' by.
    var value = this['state']['inputValue'].toLowerCase();

    // The subset of permissions that have valid roles.
    var permsWithValidRoles = [];

    // If the user is an editor or admin on this project then show the entry.
    for(var index in permissions){

      var role = permissions[index]['role'];
      if((role == 'editor') || (role == 'administrator')){

        permsWithValidRoles.push(permissions[index]);
      }
    }

    /*
     * Higher up in the render chain there is a check for permissions with valid 
     * roles (editor or admininstrator). You should not be able to add a new 
     * file with out the proper permissions (The UI should be disabled). So, 
     * 'permsWithVaildRoles' should have at least one entry if you are here.
     */

    // The subset of permissions that match the search.
    var permsMatchingSearch = [];

    // If the search value is longer than 3 characters, then filter permissions.
    if(value['length'] >= 3){

      for(var index in permsWithValidRoles){

        var projectName = permsWithValidRoles[index]['projectName'];
        if(projectName.indexOf(value) !== -1){

          permsMatchingSearch.push(permissions[index]);
        }
      }
    }
    else{

      permsMatchingSearch = permsWithValidRoles;
    }

    // In the event of no matches, say so...
    if(permsMatchingSearch == 0){

      permsMatchingSearch = [{

        'projectId': null,
        'role': null,
        'projectName': 'no matches...'
      }]
    }

    return permsMatchingSearch;
  }

  renderProjectTray(){

    var permissions = this.filterPermissions();
    if(this['state']['projectTrayActive']){

      return (
        
        <div className='project-dropdown-tray'>

          { permissions.map((permission, index)=>{

            var trayEntryProps = {

              'className': 'project-tray-entry',
              'key': index,
              'onClick': this.trayEntrySelected.bind(this),
              'data-name': permission['projectName'],
              'data-id': permission['projectId'],
              'data-role': permission['role']
            };

            return ( 

              <div { ...trayEntryProps }>

                { permission['projectName'] }
              </div> 
            );
          }) }
        </div>
      );
    }
    else{

      return '';
    }
  }

  render(){

    var listEntryProjectGroup = {

      'className': 'list-entry-project-group',
      'onMouseLeave': this['disableTray'].bind(this)
    };

    var projectSearchGroup = {

      'className': 'project-search-group',
      'title': 'The project this file belongs to.'
    };

    var permNameInputProps = {

      'className': 'project-search-input',
      'value': this['state']['inputValue'],
      'onChange': this['updateFileName'].bind(this)
    };

    var projectDropdownButton = {

      'className': 'project-dropdown-button',
      'onClick': this['toggleTray'].bind(this)
    };

    var listEntryStatus = {

      'className': 'list-entry-status project-permission-field',
      'title': 'Your project permission for this file.'
    };

    return (

      <td { ...listEntryProjectGroup }>

        <div { ...projectSearchGroup }>

          <input { ...permNameInputProps } />
          <button { ...projectDropdownButton }>

            <span className='glyphicon glyphicon-triangle-bottom'></span>
          </button>
          { this.renderProjectTray() }
        </div>
        <div { ...listEntryStatus }>

          <span className='light-text'>

            { this['state']['projectRole'] }
          </span>
        </div>
      </td>
    );
  }
}