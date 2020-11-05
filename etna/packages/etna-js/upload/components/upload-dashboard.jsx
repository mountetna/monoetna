import {selectUploads} from '../../selectors/directory-selector';

const COLUMNS = [
  {name: 'name', width: '60%'},
  {name: 'status', width: '90px', hide: true},
  {name: 'updated', width: '30%'},
  {name: 'size', width: '10%'},
  {name: 'control', width: '100px', hide: true}
];

const COLUMN_WIDTHS = COLUMNS.reduce((widths, column) => {
  widths[column.name] = column.width;
  return widths;
}, {});
