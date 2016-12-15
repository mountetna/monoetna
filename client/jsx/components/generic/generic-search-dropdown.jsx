import * as React from 'react';

export default class GenericSearchDropdown extends React.Component{

  constructor(props){

    super(props);

    this['state'] = {

      'inputValue': '',
      'trayActive': false,
      'searchEnabled': false
    };
  }

  setInputClass(){

    var className = 'admin-entry-input';
    if(this['props']['editActive']){

      className = 'admin-entry-input';

      if(this['state']['trayActive']){

        className = 'admin-entry-input admin-entry-input-active';
      }
      else{

        className = 'admin-entry-input';
      }
    }
    else{

      className = 'admin-entry-input-inactive';
    }

    return className;
  }

  setInputDisabled(){

    var disabled = false
    if(this['props']['editActive']){

      disabled = false;

      if(this['state']['searchEnabled']){

        disabled = false;        
      }
      else{

        disabled = true;
      }
    }
    else{

      disabled = true;
    }

    return disabled;
  }

  setButtonStyle(){

    var btnStyle = { 'display': 'block' };
    if(this['props']['editActive']){

      btnStyle = { 'display': 'block' };
    }
    else{

      btnStyle = { 'display': 'none' };
    }

    return btnStyle;
  }

  setDropdownStyle(){

    var dropdownStyle = { 'display': 'block' };
    if(this['props']['editActive']){

      dropdownStyle = { 'display': 'block', 'zIndex': 10000 };

      if(this['state']['trayActive']){

        dropdownStyle = { 'display': 'block', 'zIndex': 10000 };
      }
      else{

        dropdownStyle = { 'display': 'none', 'zIndex': 0 };
      }
    }
    else{

      dropdownStyle = { 'display': 'none', 'zIndex': 0 };
    }

    return dropdownStyle;
  }

  entrySelectedByClick(event){

    this.setState({ 

      'inputValue': event['target'].getAttribute('data-val'),
      'trayActive': false
    });
  }

  toggleDropdown(event){

    var trayActive = (this['state']['trayActive']) ? false : true;
    this.setState({ 'trayActive': trayActive });
  }

  openDropdownFromClick(event){

    if(!this['state']['searchEnabled']){

      this.toggleDropdown(event);
    }
  }

  updateInputValue(event){

    var value = event['target']['value'];
    if(this['state']['searchEnabled'] && value['length'] >= 3){

      var state = { 'inputValue': value, 'trayActive': true };
    }
    else if(this['state']['searchEnabled'] && value['length'] < 3){

      var state = { 'inputValue': value, 'trayActive': false };
    }
    else{

      var state = { 'inputValue': value };
    }

    this.setState(state);
  }

  // We do a simple string pattern match and return an entry if appropriate.
  matchAndAdd(value, entry, index, id){

    var text = entry;
    var entry = String(entry).toLowerCase();
    if(entry.indexOf(value) !== -1){

      var entryProps = {

        'className': 'search-dropdown-tray-entry',
        'key': (entry +'-'+ index),
        'data-id': id,
        'data-val': entry,
        'onClick': this['entrySelectedByClick'].bind(this)
      };

      return (

        <button { ...entryProps }> 

          { text }
        </button>
      );
    }
    else{

      return null;
    }
  }

  addEmpty(){

    return (

      <div className='search-dropdown-tray-empty'> 

        <i>

          { 'no matching entries...' }
        </i>
      </div>
    );
  }

  addSuggestion(){

    return (

      <div className='search-dropdown-tray-empty'> 

        <i>

          { 'type to search...' }
        </i>
      </div>
    );
  }

  // This function is overwritten by the inheriting class.
  render(){

    return (<div></div>);
  }
}