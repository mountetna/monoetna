import * as React from 'react';

import UploadMeter from './upload-meter';
import UploadControl from './upload-control';

const UploadEntryColumn = ({className, widths, children}) => (
  <div
    className={`list-entry-column-group ${className}`}
    style={{flexBasis: widths[className]}}
  >
    {children}
  </div>
);

export default UploadEntry = ({upload, widths}) => (
  <div className='list-entry-group upload'>
    <UploadEntryColumn className='name' widths={widths}>
      <span className='list-entry-file-name'>{upload.file_name}</span>
    </UploadEntryColumn>
    <UploadEntryColumn className='status' widths={widths} />
    <UploadEntryColumn className='updated' widths={widths}>
      <UploadMeter upload={upload} />
    </UploadEntryColumn>
    <UploadEntryColumn className='size' widths={widths} />
    <UploadEntryColumn className='control' widths={widths}>
      <UploadControl upload={upload} />
    </UploadEntryColumn>
  </div>
);
