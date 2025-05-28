import React, { useContext, useEffect } from 'react';
import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
import TextField from '@material-ui/core/TextField';
import Select from '@material-ui/core/Select';
import MenuItem from '@material-ui/core/MenuItem';
import Switch from '@material-ui/core/Switch';

import {PickFolder} from 'etna-js/components/metis_exploration';

import {makeStyles, Theme} from '@material-ui/core/styles';

import {ScriptItem} from '../polyphemus';
import SelectAttribute from '../select-attribute';
import TestFileMatch from './test-file-match';
import SelectAttributes from '../select-attributes';
import {MagmaContext} from 'etna-js/contexts/magma-context';

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
  <MenuItem value=''><em>{'Table format'}</em></MenuItem>
  <MenuItem value='auto-detect'>auto-detect</MenuItem>
  <MenuItem value='tsv'>tsv</MenuItem>
  <MenuItem value='csv'>csv</MenuItem>
</Select>;

export const BlankTable = ({value, update, modelName, classes}:ScriptItem) => {
  const {models} = useContext(MagmaContext);
  const parentAttribute: string = models ? models[modelName]?.template?.parent : '__error__'
  const isTable: boolean = parentAttribute!='project' && models[parentAttribute]?.template?.attributes[modelName].attribute_type=='table'

  if (value==undefined) {
    update(true)
  }

  const the_switch = <Grid item>
      <Switch
        checked={isTable && value} onChange={()=>update(!value)}
        disabled={!isTable}
        color="primary"
        />
    </Grid>

  return <Typography component="div">
      {
        isTable ?
          <Grid item container direction='row'
            component="label"
            alignItems="center"
            spacing={1}>
            <Grid item>
              Append to Previous Values
            </Grid>
            <Grid item>
              {the_switch}
            </Grid>
            <Grid item>
              Clear Previous Values
            </Grid>
          </Grid> :
          <Grid item container direction='row' alignItems="center">
            <Grid item>
              {the_switch}
            </Grid>
            <Grid item>
              Table Models Only
            </Grid>
          </Grid>
      }
  </Typography>
}

export const ValuesToIgnore = ({value, update, modelName, classes}:ScriptItem) => {
  function displayTargetValues(value: string) {
    // wrap each in 's unless empty string, and stress empty string
    const targets = value.split(',')
      .map(function(x) {
        return ("'" + x + "'").replace(/^''$/, 'empty string')
      })
    const start = `${targets.length} Target${(targets.length > 1) ? 's' : ''}`
    return `${start}: ` + (targets.join(', '))
  }
  return <Grid item container>
    <Grid item xs={6}>
      <TextField
        placeholder={'comma separated list with space respected'}
        InputLabelProps={{shrink: true}}
        fullWidth
        value={value}
        helperText={displayTargetValues(value as string)}
        onChange={
          (e: React.ChangeEvent<HTMLInputElement>) => update(e.target.value)
        }
      />
    </Grid>
  </Grid>;
}

export const ColumnMap = ({value, update, modelName, classes}:ScriptItem) => {

  const {models} = useContext(MagmaContext);
  const model = models ? models[modelName] : null;

  const idAttribute: string = models?.[modelName]?.template?.identifier || '__error__'
  const parentAttribute: string = models?.[modelName]?.template?.parent || '__error__'
  // Determine if table.
  const isTable: boolean = parentAttribute!='project' && models[parentAttribute]?.template?.attributes[modelName].attribute_type=='table'
  const autoAttribute: string = isTable ? parentAttribute : idAttribute

  useEffect( () => {
    if (autoAttribute && !(autoAttribute in value)) {
      update({
        [autoAttribute]: autoAttribute,
        ...value
      });
    }
  }, []);

  const attributesChosen = Object.keys(value).filter((val)=>val!=autoAttribute)
  const attributesToFilter = [autoAttribute]

  function onAttributeSelection(userChoice: string[]) {
    const nextValues = {...value}
    // New attributes
    const newKeys = userChoice.filter( v => !Object.keys(value).includes(v))
    if (newKeys.length > 0) {
      for (let add of newKeys) {
        nextValues[add]= add
      }
    }
    // Removed attributes
    const removedKeys = Object.keys(value).filter( v => ![...userChoice, autoAttribute].includes(v))
    if (removedKeys.length > 0) {
      for (let remove of removedKeys) {
        delete nextValues[remove]
      }
    }
    update(nextValues)
  }

  const columnNameEditor = (attributeName: string) => {
    return <TextField
      fullWidth
      value={value[attributeName]}
      onChange={(event: any) => {
        const nextValues = {...value}
        nextValues[attributeName] = event.target.value
        update(nextValues)
      }}
    />
  }

  return <Grid item container direction='column'>
    {/* Attributes Selection */}
    <Grid item>
      <SelectAttributes
        modelName={modelName}
        value={attributesChosen}
        update={onAttributeSelection}
        reportError={attributesChosen.length<1}
        placeholder='At least one additional attribute is required'
        // don't show identifier/parent & disallow file and file_collection type attributes
        filter = {(a) => !['file', 'file_collection'].includes(a.attribute_type) && !attributesToFilter.includes(a.name)}
      />
    </Grid>
    {/* Column Name Editor Table */}
    <Grid item container>
      <Grid key='titles' item container direction='row' style={{paddingTop:6, paddingBottom:4}}>
        <Grid item xs={4}><em>Attribute</em></Grid>
        <Grid item xs={8}><em>Column name in data_frame files</em></Grid>
      </Grid>
      {
        Object.keys(value).map(
          (attribute_name) => {
            let attribute_label = attribute_name
            if (attribute_name==idAttribute) attribute_label += ' (Identifier)'
            if (attribute_name==parentAttribute) attribute_label += ' (Parent)'
            return <Grid
              key={attribute_name} item container
              direction='row'
              >
              <Grid item xs={4}
                className={attribute_name==autoAttribute ? undefined : classes.strikeout}
                onClick={attribute_name==autoAttribute ? undefined : () => {
                  let newValues = { ...value };
                  delete newValues[attribute_name];
                  update(newValues);
                }
              }>{attribute_label}</Grid>
              <Grid item xs={8}>{columnNameEditor(attribute_name)}</Grid>
            </Grid>
          }
        )
      }
    </Grid>
  </Grid>;
}

// export const ExtractedColumns = ({value,update}:ScriptItem) => <TextField
//   placeholder='Comma-separated list'
//   fullWidth
//   value={Array.isArray(value) ? value.join(', ') : value}
//   onChange={
//     (e: React.ChangeEvent<HTMLInputElement>) => update(e.target.value.split(/,\s*/))
//   }
// />;
