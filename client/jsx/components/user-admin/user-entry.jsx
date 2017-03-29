import * as React from 'react';
import GenericAdminEntry from './generic-admin-entry';

export default class UserEntry extends GenericAdminEntry{

  render(){

    var adminEntryProps = {

      'className': 'admin-edit-entry-group',
      'onMouseEnter': this['showControlGroup'].bind(this),
      'onMouseLeave': this['hideControlGroup'].bind(this)
    };

    var emailEntryProps = {

      'className': 'admin-entry-input',
      'value': this['props']['user']['email'],
      'title': this['props']['user']['email']
    };

    var firstNameEntryProps = {

      'className': 'admin-entry-input',
      'value': this['props']['user']['firstName'],
      'title': this['props']['user']['firstName']
    };

    var lastNameEntryProps = {

      'className': 'admin-entry-input',
      'value': this['props']['user']['lastName'],
      'title': this['props']['user']['lastName']
    };

    if(!this['state']['editActive']){

      emailEntryProps['className'] = 'admin-entry-input-inactive';
      emailEntryProps['disabled'] = true;

      firstNameEntryProps['className'] = 'admin-entry-input-inactive';
      firstNameEntryProps['disabled'] = true;

      lastNameEntryProps['className'] = 'admin-entry-input-inactive';
      lastNameEntryProps['disabled'] = true;
    }

    return (

      <tr { ...adminEntryProps }>

        <td>

          <input { ...emailEntryProps } />
        </td>
        <td>

          <input { ...firstNameEntryProps } />
        </td>
        <td>

          <input { ...lastNameEntryProps } />
        </td>
        <td>
  
          { this.renderEditControlGroup() }
        </td>
      </tr>
    );
  }
}