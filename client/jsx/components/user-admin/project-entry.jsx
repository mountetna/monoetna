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

      'className': 'admin-entry-input',
      'value': this['props']['project']['projectName'],
      'title': this['props']['project']['projectName']
    };

    var groupNameProps = {

      'className': 'admin-entry-input',
      'value': this['props']['project']['groupName'],
      'title': this['props']['project']['groupName']
    }

    if(!this['state']['editActive']){

      projectNameProps['className'] = 'admin-entry-input-inactive';
      projectNameProps['disabled'] = true;

      groupNameProps['className'] = 'admin-entry-input-inactive';
      groupNameProps['disabled'] = true;
    }

    return (

      <tr { ...adminEntryProps }>

        <td>

          <input { ...projectNameProps } />
        </td>
        <td>

          <input { ...groupNameProps } />
        </td>
        <td>

          { this.renderEditControlGroup() }
        </td>
      </tr>
    );
  }
}