import React, { useState, useEffect, useCallback } from 'react';
import {makeStyles} from '@material-ui/core/styles';
import CardActions from '@material-ui/core/CardActions';
import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
import TextField from '@material-ui/core/TextField';
import Button from '@material-ui/core/Button';
import HistoryIcon from '@material-ui/icons/HistoryRounded';
import Tooltip from '@material-ui/core/Tooltip';
import IconButton from '@material-ui/core/IconButton';

import ConfigScript from './config-script';
import EtlPane, { EtlPaneHeader } from './etl-pane';
import RevisionHistory from './revision-history';
import { formatTime } from './run-state';

const useStyles = makeStyles( theme => ({
  savebar: {
    justifyContent: 'space-between',
    width: 'auto'
  },
  revisions: {
    flex: '1 1 auto',
    width: 'auto',
    justifyContent: 'flex-end'
  }
}));

const ConfigurePane = ({name, project_name, selected, config, schema, update}:{name:string, project_name:string, selected:string|null, config:any, schema:any,update:Function}) => {
  const classes = useStyles();

  let [ origScript, setOrigScript ] = useState('');
  let [ editedScript, setEditedScript ] = useState('');
  let [ comment, setComment ] = useState('');
  let [ showRevisions, setShowRevisions ] = useState(false);

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
        <Grid container spacing={1} item className={classes.revisions}>
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
    <ConfigScript script={editedScript} schema={schema} onChange={setEditedScript} />
  </EtlPane>
}

export default ConfigurePane;
