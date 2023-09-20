import React from 'react';
import Grid from '@material-ui/core/Grid';
import TextField from '@material-ui/core/TextField';
import Select from '@material-ui/core/Select';
import MenuItem from '@material-ui/core/MenuItem';

import AddDialog from 'etna-js/components/add-dialog';
import {PickFolder} from 'etna-js/components/metis_exploration';

import {makeStyles, Theme} from '@material-ui/core/styles';

import {ScriptItem} from '../polyphemus';
import SelectAttribute from '../select-attribute';
import TestFileMatch from './test-file-match';

export const useMetisFormStyles = makeStyles((theme: Theme) => ({
  form: {
    height: 'calc(100vh - 375px)',
    flexWrap: 'nowrap'
  },
  tabs: {
    flex: '1 1 auto',
    maxHeight: '50px',
    border: '1px solid #ccc',
    borderBottom: 'none'
  },
  tab_scroll: {
    flexGrow: 0,
    maxWidth: 'calc(100% - 120px)',
    flexBasis: 'calc(100% - 120px)'
  },
  tab_buttons: {
    flex: '1 1 auto',
    width: 'auto'
  },
  tab_pane: {
    flex: '1 1 auto',
    border: '1px solid #ccc',
    overflow: 'hidden',
    overflowY: 'auto'
  },
  bucket_name: {
    padding: '10px',
    width: '500px'
  },
  model: {
    padding: '5px',
    paddingRight: '15px'
  },
  model_row: {
    minHeight: '47px',
    '&:not(:last-of-type)': {
      borderBottom: '1px solid #eee',
      marginBottom: '5px'
    }
  },
  model_row_name: {
    marginTop: '7px'
  },
  page_size: {
    color: '#666',
    padding: '0px 15px'
  },
  script: {
    background: '#eee',
    marginTop: '15px',
    border: '1px solid #ccc',
    borderRadius: '2px',
    position: 'relative'
  },
  script_header: {
    borderBottom: '1px solid #ccc',
    background: '#ccc'
  },
  number: {
    color: '#888',
    position: 'absolute',
    transform: 'translate(-100%, 0)',
    left: '-10px',
    top: '10px'
  },
  attribute_name: {
    overflow: 'hidden',
    textOverflow: 'ellipsis',
    padding: '0px 5px'
  },
  script_items: {
    padding: '5px'
  },
  strikeout: {
    '&:hover': {
      textDecoration: 'line-through',
      cursor: 'pointer',
      background: '#eee'
    }
  },
  test: {
    padding: '0px 10px'
  }
}));

export const DefaultItem = ({value, update}:ScriptItem) => <TextField
  fullWidth
  value={value}
  onChange={ e => update(e.target.value as string) }
/>;

export const FileMatch = ({value,update,projectName,bucketName,script,classes}:ScriptItem) =>  {
  return <Grid item container>
    <Grid item xs={3}>
      <TextField
        placeholder='Glob matching file, e.g. *.fastq.gz'
        fullWidth
        value={value}
        onChange={
          (e: React.ChangeEvent<HTMLInputElement>) => update(e.target.value)
        }
      />
    </Grid>
    <TestFileMatch className={classes.test} projectName={projectName} bucketName={bucketName} script={script}/>
  </Grid>;
}

export const FolderPath = ({value,update,projectName,bucketName}:ScriptItem) => <Grid item xs={3}>
  <PickFolder
    basePath=''
    label={undefined}
    path={value}
    setPath={update}
    project_name={projectName}
    bucket={bucketName}
  />
</Grid>;

export const AttributeName = ({type,value,update,modelName}:ScriptItem) => <SelectAttribute
  value={value}
  update={update}
  modelName={modelName}
  filter={ (a:any) => a.attribute_type == type }
/>;


export const TableFormat = ({value,update}:ScriptItem) => <Select displayEmpty value={value} onChange={e => update(e.target.value)}>
  <MenuItem value=''><em>Table format</em></MenuItem>
  <MenuItem value='tsv'>tsv</MenuItem>
  <MenuItem value='csv'>csv</MenuItem>
</Select>;

export const ColumnMap = ({value,update,classes}:ScriptItem) => <Grid>
  <AddDialog
    title='Map a column'
    content={
      <>Map a column name in the data_frame to an attribute_name</>
    }
    buttonText='ADD MAPPING'
    update={ (column_name:string, mapped_column_name:string) => update(
      { ...value, [column_name]: mapped_column_name }
    )}
    mask={ v => v }
    mask2={ v => v }
    placeholders={[ 'Column Name', 'Mapped Column Name' ]}
  />
  <Grid item container>
    {
      Object.keys(value).map(
        column_name => <Grid
          key={column_name} item container
          className={classes.strikeout}
          onClick={
            () => {
              let newValue = { ...value };
              delete newValue[column_name];
              update(newValue);
            }
          }>
          <Grid item xs={2}>{column_name}</Grid>
          <Grid item xs={8}>{value[column_name]}</Grid>
        </Grid>
      )
    }
  </Grid>
</Grid>;

export const ExtractedColumns = ({value,update}:ScriptItem) => <TextField
  placeholder='Comma-separated list'
  fullWidth
  value={Array.isArray(value) ? value.join(', ') : value}
  onChange={
    (e: React.ChangeEvent<HTMLInputElement>) => update(e.target.value.split(/,\s*/))
  }
/>;