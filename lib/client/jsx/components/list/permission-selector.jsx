import * as React from 'react';
import GenericSearchDropdown from '../generic/generic-search-dropdown';

export default class PermissionSelector extends GenericSearchDropdown{
  componentDidMount(){
    this.setState({
      searchEnabled: true
    });
  }

  addEntries(){
    // The value of the input to 'search' by.
    let value = this.state.inputValue;

    /*
     * If there are more than three characters in the input, then 'search'.
     * This prevents unnecessary processing for the UI.
     */
    if (value.length >= 3){
      // Generate a list of  entries to display in the dropdown.
      let entries = [];

      // The list of user permissions to be 'searched' over.
      let { permissions } = this.props;

      permissions.forEach((permission,i) => {
        let { project_name, role } = permission;
        let entry = this.matchAndAdd(value, project_name, i, entries.length);
        
        if (entry != null && role != 'viewer') entries.push(entry);
      })

      // If there are no matching entries then add a blank entry that says so.
      if(entries.length == 0){
        entries.push(this.addEmpty());
      }

      // Return a mapping to render the entries.
      return entries;
    }
    else{
      return this.addSuggestion();
    }
  }

  entrySelected(){
    let { inputValue }  = this.state;
    let { permissions } = this.props;
    let permission = permissions.find( perm => perm.project_name == inputValue );

    if(permission != null){
      this.props.callbacks.projectSelected(permission);
    }
  }

  setRole(){
    let { inputValue } = this.state;
    let { permissions } = this.props;
    let permission = permissions.find(perm => perm.project_name == inputValue)

    return permission ? permission.role : '';
  }

  render(){
    let listEntryProjectGroup = {
      className: 'list-entry-project-group'
    };

    let dropdownGroupProps = {
      className: 'list-project-search-dd-group',
      title: 'The project this file belongs to.'
    };

    let dropdownInputProps = {
      className: this.setInputClass(),
      disabled: this.setInputDisabled(),
      value: this.state.inputValue,
      onChange: this.updateInputValue.bind(this),
      onKeyUp: this.selectByKeyboard.bind(this),
      ref: (component)=>{ this.dropdownInput = component }
    };

    let dropdownBtnProps = {
      className: 'list-project-search-dd-btn',
      onClick: this.toggleDropdown.bind(this),
      style: this.setButtonStyle()
    };

    let dropdownTrayProps = {
      className: 'list-project-search-dd-tray',
      style: this.setDropdownStyle(),
      ref: (component)=>{ this.dropdownTrayComponent = component }
    };

    let listEntryStatus = {
      className: 'list-entry-status list-project-permission-field',
      title: 'Your project permission for this file.'
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
