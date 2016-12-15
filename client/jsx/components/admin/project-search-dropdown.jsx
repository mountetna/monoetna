import * as React from 'react';
import GenericSearchDropdown from '../generic/generic-search-dropdown';

export default class ProjectSearchDropdown extends GenericSearchDropdown{

  componentDidMount(){

    this.setState({ 

      'searchEnabled': true,
      'inputValue': this['props']['permission']['projectName'],
    });
  }

  addEntries(){

    // The list of entries to be 'searched' over.
    var projects = this['props']['projects'];

    // The value of the input to 'search' by.
    var value = this['state']['inputValue'].toLowerCase();

    /*
     * If there are more than three characters in the input, then 'search'.
     * This prevents unnecessary processing for the UI.
     */
    if(value['length'] >= 3){

      // Generate a list of entries to display in the dropdown.
      var entries = [];
      for(var a = 0; a < projects['length']; ++a){

        var name = projects[a]['projectName'];
        var id = projects[a]['id'];
        var entry = this.matchAndAdd(value, name, a, id);
        if(entry != null) entries.push(entry);
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

  render(){

    // Return the input value to the parent.
    var inputVal = this['state']['inputValue'];
    this['props']['callbacks'].updateFromInput('projectName', inputVal);

    var dropdownGroupProps = {

      'className': 'perm-project-search-dd-group'
    };

    var dropdownInputProps = {

      'className': this.setInputClass(),
      'disabled': this.setInputDisabled(),
      'value': this['state']['inputValue'],
      'onChange': this['updateInputValue'].bind(this)
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

    return (

      <div { ...dropdownGroupProps }>
        
        <input { ...dropdownInputProps } />
        <button { ...dropdownBtnProps } >

          <span className='glyphicon glyphicon-triangle-bottom'></span>
        </button>

        <div { ...dropdownTrayProps }>

          { this.addEntries() }
        </div>
      </div>
    );
  }
}