import React, { useState, useCallback } from 'react';
import EtlPane, {EtlPaneHeader} from './etl-pane';
import {makeStyles} from '@material-ui/core/styles';
import { Controlled } from 'react-codemirror2';
import Typography from '@material-ui/core/Typography';
import MenuItem from '@material-ui/core/MenuItem';
import Select from '@material-ui/core/Select';
import TextField from '@material-ui/core/TextField';
import Grid from '@material-ui/core/Grid';
import Button from '@material-ui/core/Button';

const useStyles = makeStyles( theme => ({
  table: {
  },
  tablerow: {
    '&:not(:last-of-type)': {
      borderBottom: '1px solid #ccc'
    },
    padding: '5px',
    fontSize: '0.8rem'
  }
}));

const SecretsPane = ({selected, keys, update, secrets}:{
  selected:string|null,
  keys:string[],
  update:Function,
  secrets:any
}) => {
  const classes = useStyles();

  const [ secret, setSecret ] = useState('');
  const [ value, setValue ] = useState('');

  return <EtlPane mode='secrets' selected={selected}>
    <EtlPaneHeader title='Secrets'/>
    <Grid container className={classes.table}>
      {
        keys.map( key => <Grid item className={classes.tablerow} container key={key}>
          <Grid item xs={4}>{key}</Grid>
          <Grid item xs={4}>{secrets[key] ? '●●●' : <em>none</em>}</Grid>
        </Grid>)
      }
    </Grid>
    <Grid container spacing={1}>
      <Grid item>
        <Select displayEmpty value={secret} onChange={e => setSecret(e.target.value as string)}>
          <MenuItem disabled value=''>Set secret</MenuItem>
          {
            keys && keys.map( key => <MenuItem key={key} value={key}>{key}</MenuItem>)
          }
        </Select>
      </Grid>
      {
        secret && <React.Fragment>
          <Grid item>
            <TextField onChange={ e => setValue(e.target.value) } fullWidth value={value} placeholder='New value'/>
          </Grid>
          <Grid item>
             <Button onClick={ () => update({secrets: { [secret]: value }}) }>Save</Button>
          </Grid>
        </React.Fragment>
      }
    </Grid>
  </EtlPane>
}

export default SecretsPane;
