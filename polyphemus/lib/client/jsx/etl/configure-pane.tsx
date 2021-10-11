import React, { useState, useEffect, useCallback } from 'react';
import {makeStyles} from '@material-ui/core/styles';
import CardActions from '@material-ui/core/CardActions';
import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
import TextField from '@material-ui/core/TextField';
import Button from '@material-ui/core/Button';
import HistoryIcon from '@material-ui/icons/HistoryRounded';
import CodeIcon from '@material-ui/icons/Code';
import CodeOffIcon from '@material-ui/icons/CodeOff';
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

const FORMS = {
  redcap: RedcapForm
}

const ConfigurePane = ({name, project_name, selected, config, job, update}:{name:string, project_name:string, selected:string|null, config:any, job:Job|null,update:Function}) => {
  const classes = useStyles();

  let [ origScript, setOrigScript ] = useState('');
  let [ editedScript, setEditedScript ] = useState('');
  let [ comment, setComment ] = useState('');
  let [ showRevisions, setShowRevisions ] = useState(false);
  let [ showJson, setShowJson ] = useState(false);

  let JobForm = job ? FORMS[job.name] : null;

  useEffect( () => {
    origScript = JSON.stringify(config,null,2);
    setOrigScript(origScript);
    setEditedScript(origScript);
  }, [config]);

  return <EtlPane mode='configure' selected={selected}>
    <EtlPaneHeader title='Configuration'>
      {
        origScript != editedScript && 
          <Grid spacing={1} container className={classes.savebar}>
            <Grid item><TextField style={{width:300}} value={comment} onChange={ e => setComment(e.target.value) } placeholder='Revision comment'/></Grid>
            <Grid item><Button disabled={ comment == '' } onClick={ () => {
              update({ comment, config: JSON.parse(editedScript) });
              setComment('');
            } } >Save</Button></Grid>
            <Grid item><Button onClick={() => setEditedScript(origScript)} color='secondary'>Reset</Button></Grid>
          </Grid>
      }
      <Grid container spacing={1} item className={classes.buttons}>
        {
          JobForm && <Tooltip title={ showJson ? 'show form' : 'show json' }>
            <IconButton onClick={() => setShowJson(!showJson) } size='small' aria-label='revision history'
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
            config => {
              setEditedScript(JSON.stringify(config,null,2));
              setShowRevisions(false)
            }
          }
          onClose={() => setShowRevisions(false)}
        />
      </Grid>
    </EtlPaneHeader>
    {
      (!JobForm || showJson)
        ? <ConfigScript script={editedScript} schema={job ? job.schema : null} onChange={setEditedScript} />
        : <JobForm config={config} />
    }
  </EtlPane>
}

export default ConfigurePane;
