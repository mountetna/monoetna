import * as React from 'react';
import GenericSearchDropdown from '../generic/generic-search-dropdown';

export default class PermissionSelector extends GenericSearchDropdown{

  componentDidMount(){

    this.setState({

      'searchEnabled': true
    });
  }

  addEntries(){

    // The value of the input to 'search' by.
    var value = this['state']['inputValue'].toLowerCase();

    /*
     * If there are more than three characters in the input, then 'search'.
     * This prevents unnecessary processing for the UI.
     */
    if(value['length'] >= 3){

      // Generate a list of  entries to display in the dropdown.
      var entries = [];

      // The list of user permissions to be 'searched' over.
      var permissions = this['props']['permissions'];

      for(var a = 0; a < permissions['length']; ++a){

        var name = permissions[a]['projectName'];
        var id = permissions[a]['id'];
        var entry = this.matchAndAdd(value, name, a, id, entries['length']);
        
        if(entry != null && permissions[a]['role'] != 'viewer'){

          entries.push(entry);
        }
      }

      // If there are no matching entries then add a blank entry that says so.
      if(entries.length == 0){

        entries.push(this.addEmpty());
      }

      // Return a mapping to render the entries.
      return entries.map((entry, index)=>{

        return entry;
      });
    }
    else{

      return this.addSuggestion();
    }
  }

  entrySelected(){

    var projectName = this['state']['inputValue'];
    var permissions = this['props']['permissions'];
    var permission = null;

    for(var a = 0; a < permissions['length']; ++a){

      if(projectName == permissions[a]['projectName'].toLowerCase()){

        permission = permissions[a];
      }
    }

    if(permission != null){

      this['props']['callbacks'].projectSelected(permission);
    }
  }

  setRole(){

    var projectName = this['state']['inputValue'];
    var permissions = this['props']['permissions'];
    var role = '';
    for(var a = 0; a < permissions['length']; ++a){

      if(permissions[a]['projectName'] == projectName){

        role = permissions[a]['role'];
        break;
      }
    }
    return role;
  }

  render(){

    var listEntryProjectGroup = {

      'className': 'list-entry-project-group'
    };

    var dropdownGroupProps = {

      'className': 'list-project-search-dd-group',
      'title': 'The project this file belongs to.'
    };

    var dropdownInputProps = {

      'className': this.setInputClass(),
      'disabled': this.setInputDisabled(),
      'value': this['state']['inputValue'],
      'onChange': this['updateInputValue'].bind(this),
      'onKeyUp': this['selectByKeyboard'].bind(this),
      'ref': (component)=>{ this['dropdownInput'] = component }
    };

    var dropdownBtnProps = {

      'className': 'list-project-search-dd-btn',
      'onClick': this['toggleDropdown'].bind(this),
      'style': this.setButtonStyle()
    };

    var dropdownTrayProps = {

      'className': 'list-project-search-dd-tray',
      'style': this.setDropdownStyle(),
      'ref': (component)=>{ this['dropdownTrayComponent'] = component }
    };

    var listEntryStatus = {

      'className': 'list-entry-status list-project-permission-field',
      'title': 'Your project permission for this file.'
    };

    return (

      <td { ...listEntryProjectGroup }>

        <div { ...dropdownGroupProps }>

          <input { ...dropdownInputProps }/>
          <button { ...dropdownBtnProps }>

            <span className='glyphicon glyphicon-triangle-bottom'></span>
          </button>
          <div { ...dropdownTrayProps }>

            { this.addEntries() }
          </div>
        </div>
        <div { ...listEntryStatus }>

          <span className='light-text'>

            { this.setRole() }
          </span>
        </div>
      </td>
    );
  }
}