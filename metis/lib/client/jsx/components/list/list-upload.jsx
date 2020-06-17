import * as React from 'react';
import { ListEntryColumn, ListEntryTypeColumn, ListEntryUpdatedColumn } from './list-entry';
import UploadMeter from './upload-meter';
import UploadControl from './upload-control';

const ListUpload = ({ upload, widths }) => (
  <div className='list-entry-group upload'>
    <ListEntryTypeColumn icon='upload' widths={ widths } />
    <ListEntryColumn className='name' widths={ widths }>
      <span className='list-entry-file-name'>
        {upload.file_name}
      </span>
    </ListEntryColumn>
    <ListEntryColumn className='status' widths={widths}/>
    <ListEntryColumn className='updated' widths={widths}>
      <UploadMeter upload={ upload }/>
    </ListEntryColumn>
    <ListEntryColumn className='size' widths={widths}/>
    <ListEntryColumn className='control' widths={widths}>
      <UploadControl upload={ upload } />
    </ListEntryColumn>
  </div>
);

export default ListUpload;
