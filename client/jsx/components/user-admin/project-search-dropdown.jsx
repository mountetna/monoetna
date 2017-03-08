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

    // The value of the input to 'search' by.
    var value = this['state']['inputValue'].toLowerCase();

    /*
     * If there are more than three characters in the input, then 'search'.
     * This prevents unnecessary processing for the UI.
     */
    if(value['length'] >= 3){

      // The list of entries to be 'searched' over.
      var projects = this['props']['projects'];

      // Generate a list of  entries to display in the dropdown.
      var entries = [];
      for(var a = 0; a < projects['length']; ++a){

        var name = projects[a]['projectName'];
        var id = projects[a]['id'];
        var entry = this.matchAndAdd(value, name, a, id, entries['length']);
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

  entrySelected(){

    this['props']['callbacks'].projectSelected(this['state']['inputValue'])
  }

  render(){

    var dropdownGroupProps = {

      'className': 'perm-project-search-dd-group'
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

      'className': 'perm-project-search-dd-btn',
      'onClick': this['toggleDropdown'].bind(this),
      'style': this.setButtonStyle()
    };

    var dropdownTrayProps = {

      'className': 'perm-project-search-dd-tray',
      'style': this.setDropdownStyle(),
      'ref': (component)=>{ this['dropdownTrayComponent'] = component }
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