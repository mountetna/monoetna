import * as React from 'react';

import {authorFormat, dateFormat} from '../utils/format';

import Icon from './icon';

export const ListEntryColumn = ({className, widths, children}) => (
  <div
    className={`list-entry-column-group ${className}`}
    style={{flexBasis: widths[className]}}
  >
    {children}
  </div>
);

export const ListEntryTypeColumn = ({icon, widths}) => (
  <ListEntryColumn className='type' widths={widths}>
    <Icon icon={icon} />
  </ListEntryColumn>
);

export const ListEntryUpdatedColumn = ({obj, widths}) => (
  <ListEntryColumn className='updated' widths={widths}>
    <div className='list-entry-updated-name'>
      {dateFormat(obj.updated_at)} by {authorFormat(obj.author).name}
    </div>
  </ListEntryColumn>
);
