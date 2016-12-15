import * as React from 'react';
import GenericAdminEntry from './generic-admin-entry';

export default class ProjectEntry extends GenericAdminEntry{

  render(){

    var adminEntryProps = {

      'className': 'admin-edit-entry-group',
      'onMouseEnter': this['showControlGroup'].bind(this),
      'onMouseLeave': this['hideControlGroup'].bind(this)
    };

    var projectNameProps = {

      'className': 'admin-entry-input-inactive',
      'value': this['props']['project']['projectName'],
      'title': this['props']['project']['projectName']
    };

    return (

      <tr { ...adminEntryProps }>

        <td>

          <input { ...projectNameProps } />
        </td>
        <td>
  
          { 'cip' }
        </td>
        <td>
  
          { this.renderEditControlGroup() }
        </td>
      </tr>
    );
  }
}