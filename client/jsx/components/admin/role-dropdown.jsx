import * as React from 'react';
import GenericSearchDropdown from '../generic/generic-search-dropdown';

export default class RoleDropdown extends GenericSearchDropdown{

  componentDidMount(){

    this.setState({ 

      'inputValue': this['props']['permission']['role']
    });
  }

  entrySelected(){

    this['props']['callbacks'].roleSelected(this['state']['inputValue']);
  }

  render(){

    // Return the input value to the parent.
    var inputVal = this['state']['inputValue'];

    var dropdownGroupProps = {

      'className': 'perm-project-search-dd-group',
      'onClick': this['openDropdownFromClick'].bind(this)
    };

    var dropdownInputProps = {

      'className': this.setInputClass(),
      'disabled': this.setInputDisabled(),
      'value': this['state']['inputValue']
    };

    var dropdownBtnProps = {

      'className': 'perm-project-search-dd-btn',
      'onClick': this['toggleDropdown'].bind(this),
      'style': this.setButtonStyle()
    };

    var dropdownTrayProps = {

      'className': 'perm-project-search-dd-tray',
      'style': this.setDropdownStyle()
    };

    var entryProps = {

      'className': 'search-dropdown-tray-entry',
      'onClick': this['entrySelectedByClick'].bind(this),
    };

    return (

      <div { ...dropdownGroupProps }>
        
        <input { ...dropdownInputProps } />
        <button { ...dropdownBtnProps } >

          <span className='glyphicon glyphicon-triangle-bottom'></span>
        </button>

        <div { ...dropdownTrayProps }>

          <button { ...entryProps } data-val='administrator'>

            { 'administrator' }
          </button>
          <button { ...entryProps } data-val='editor'>

            { 'editor' }
          </button>
          <button { ...entryProps } data-val='viewer'>

            { 'viewer' }
          </button>
        </div>
      </div>
    );
  }
}