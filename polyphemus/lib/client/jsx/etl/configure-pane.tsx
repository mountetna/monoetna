import React, { useState, useEffect, useCallback } from 'react';
import {makeStyles} from '@material-ui/core/styles';
import CardActions from '@material-ui/core/CardActions';
import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
import TextField from '@material-ui/core/TextField';
import Button from '@material-ui/core/Button';
import HistoryIcon from '@material-ui/icons/HistoryRounded';
import CodeIcon from '@material-ui/icons/Code';
import Tooltip from '@material-ui/core/Tooltip';
import IconButton from '@material-ui/core/IconButton';
import FormGroup from '@material-ui/core/FormGroup';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import Switch from '@material-ui/core/Switch';

import ConfigScript from './config-script';
import EtlPane, { EtlPaneHeader } from './etl-pane';
import RedcapForm from './redcap-form';
import RevisionHistory from './revision-history';
import { formatTime } from './run-state';

const useStyles = makeStyles( theme => ({
  savebar: {
    justifyContent: 'space-between',
    width: 'auto'
  },
  buttons: {
    flex: '1 1 auto',
    width: 'auto',
    justifyContent: 'flex-end'
  }
}));

const FORMS:{
  [key:string]: any
} = {
  redcap: RedcapForm
};

const ConfigurePane = ({name, project_name, selected, config, job, update, error}:{
  name:string,
  project_name:string,
  selected:string|null,
  config:any,
  job:Job|undefined,
  update:Function,
  error:string
}) => {
  const classes = useStyles();

  const [ origScript, setOrigScript ] = useState('');
  const [ editedScript, setEditedScript ] = useState('');
  const [ editedConfig, setEditedConfig ] = useState({});
  const [ comment, setComment ] = useState('');
  const [ showRevisions, setShowRevisions ] = useState(false);
  const [ showJson, setShowJson ] = useState(false);

  const JobForm = job ? FORMS[job.name] : null;

  useEffect( () => {
    const script = JSON.stringify(config,null,2);
    setOrigScript(script);
    setEditedScript(script);
    setEditedConfig(config);
  }, [config]);

  const showSave = showJson ? (origScript != editedScript) : (editedConfig != config);

  return <EtlPane mode='configure' selected={selected}>
    <EtlPaneHeader title='Configuration' error={error}>
      {
        showSave && 
          <Grid spacing={1} container className={classes.savebar}>
            <Grid item><TextField style={{width:300}} value={comment} onChange={ e => setComment(e.target.value) } placeholder='Revision comment'/></Grid>
            <Grid item><Button disabled={ comment == '' } onClick={ () => {
              update({ comment, config: showJson ? JSON.parse(editedScript) : editedConfig });
              setComment('');
            } } >Save</Button></Grid >
            <Grid item><Button onClick={() => showJson ? setEditedScript(origScript) : setEditedConfig(config)} color='secondary'>Reset</Button></Grid>
          </Grid>
      }
      <Grid container spacing={1} item className={classes.buttons}>
        {
          JobForm && <Tooltip title={ showJson ? 'show form' : 'show json' }>
            <IconButton disabled={showSave} onClick={() => setShowJson(!showJson) } size='small' aria-label='revision history'
              color={ showJson ? 'primary' : 'default' }>
              <CodeIcon/>
            </IconButton>
          </Tooltip>
        }
        <Tooltip title='revision history'>
          <IconButton onClick={() => setShowRevisions(true) } size='small' aria-label='revision history'>
            <HistoryIcon/>
          </IconButton>
        </Tooltip>
        <RevisionHistory
          name={name}
          project_name={project_name}
          open={showRevisions}
          config={config}
          update={
            newConfig => {
              setEditedConfig(newConfig);
              setEditedScript(JSON.stringify(newConfig,null,2));
              setShowRevisions(false);
            }
          }
          onClose={() => setShowRevisions(false)}
        />
      </Grid>
    </EtlPaneHeader>
    {
      (!JobForm || showJson)
        ? <ConfigScript script={editedScript} schema={job ? job.schema : null} onChange={setEditedScript} />
        : <JobForm project_name={project_name} config={editedConfig} job={job} update={ setEditedConfig }/>
    }
  </EtlPane>
}

export default ConfigurePane;
