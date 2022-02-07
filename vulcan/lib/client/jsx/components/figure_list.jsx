import React, {useCallback, useContext, useState, useEffect} from 'react';
import 'regenerator-runtime/runtime';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {pushLocation} from 'etna-js/actions/location_actions';

import {VulcanContext} from '../contexts/vulcan_context';
import WorkflowCard from '../components/dashboard/card';
import {workflowName} from "../selectors/workflow_selectors";
import SelectInput from 'etna-js/components/inputs/select_input'

import Grid from '@material-ui/core/Grid';
import Menu from '@material-ui/core/Menu';
import MenuItem from '@material-ui/core/MenuItem';
import Card from '@material-ui/core/Card';
import CardHeader from '@material-ui/core/CardHeader';
import CardMedia from '@material-ui/core/CardMedia';
import CardContent from '@material-ui/core/CardContent';
import Dialog from '@material-ui/core/Dialog';
import DialogContent from '@material-ui/core/DialogContent';
import DialogActions from '@material-ui/core/DialogActions';
import DialogTitle from '@material-ui/core/DialogTitle';
import Button from '@material-ui/core/Button';
import Typography from '@material-ui/core/Typography';

import IconButton from '@material-ui/core/IconButton';
import MoreVertIcon from '@material-ui/icons/MoreVert';

import {makeStyles} from '@material-ui/core/styles';
import { json_get } from 'etna-js/utils/fetch';

const createFigureStyles = makeStyles( theme => ({
  none: {
    color: '#f44'
  },
  workflows: {
  },
  selected: {
    border: '1px solid red',
    background: 'rgba(100,0,0,0.1)',
    boxShadow: '0 0 0 15px #fee, 0 0 4px 15px #faa'
  }
}));

const CreateFigure = ({handleClose, open, project_name}) => {
  let {state} = useContext(VulcanContext);
  const {workflows} = state;

  const [ selectedWorkflow, selectWorkflow ] = useState(null);
  const classes = createFigureStyles();

  const projectWorkflows = workflows ? workflows.filter( ({projects}) => projects.includes(project_name) ) : [];

  const invoke = useActionInvoker();

  const visitNewFigure = useCallback((workflow) => {
    invoke(pushLocation(`/${project_name}/figure/new/${workflowName(workflow)}`));
  }, [invoke])

  const close = useCallback( () => {
    handleClose();
    selectWorkflow(null);
  }, [handleClose] );

  return <Dialog onClose={close} open={open}>
    <DialogTitle>Select a workflow</DialogTitle>
    <DialogContent dividers>
    {
      projectWorkflows.length ? <Grid  className={classes.workflows} container direction='column'>
        {
          projectWorkflows.map((w, ind) =>
            <WorkflowCard 
            className={ selectedWorkflow == w ? classes.selected : null }
            workflow={w} key={ind} onClick={ () => selectWorkflow( selectedWorkflow == w ? null : w) }/>
          )
        }
      </Grid> : <Grid item className={classes.none}><em>No workflows</em></Grid>
    }
    </DialogContent>
    <DialogActions>
      { selectedWorkflow && <Typography variant='body1'>{selectedWorkflow.displayName}</Typography> }
      <Button variant='contained' disabled={ !selectedWorkflow } onClick={ () => visitNewFigure(selectedWorkflow) }>Create Figure</Button>
    </DialogActions>
  </Dialog>
}

const figureStyles = makeStyles( theme => ({
  figure: {
    maxWidth: 350
  },
  image: {
    cursor: 'pointer'
  }
}));

const Figure = ({figure}) => {
  const invoke = useActionInvoker();
  let {state} = useContext(VulcanContext);
  const {workflows} = state;

  const workflow = workflows ? workflows.find( w => w.name == figure.workflow_name ) : null

  const classes = figureStyles();

  const visitFigure = useCallback((figure) => {
    invoke(pushLocation(`/${figure.project_name}/figure/${figure.id}`));
  }, [invoke]);

  const [ menuAnchor, setMenuAnchor ] = useState(null);

  const handleClose = () => { setMenuAnchor(null); }

  return <Card className={classes.figure}>
    <Menu
      id="simple-menu"
      open={Boolean(menuAnchor)}
      anchorEl={ menuAnchor }
      onClose={ handleClose }
    >
        <MenuItem onClick={handleClose}>Copy</MenuItem>
        <MenuItem onClick={handleClose}>Rename</MenuItem>
        <MenuItem onClick={handleClose}>Remove</MenuItem>
    </Menu>
    <CardHeader
      title={figure.title}
      titleTypographyProps={ { variant: 'h6' } }
      subheader={figure.workflow_name.replace('.cwl','')}
      action={
        <IconButton onClick={ (e) => setMenuAnchor(e.currentTarget) }>
          <MoreVertIcon />
        </IconButton>
      }/>
    <CardMedia
      className={classes.image}
      onClick={ () => visitFigure(figure) }
      component="img"
      height="140"
      image={`/images/${workflow ? workflow.image : 'default.png'}`}
      title={figure.title}
    />
  </Card>
}

const figureListStyles = makeStyles( theme => ({
  title: {
    padding: '10px 0px',
    color: '#444'
  },
  figures: {
    padding: '15px'
  }
}));

export default function FigureList({project_name}) {
  const invoke = useActionInvoker();
  const [ figures, setFigures ] = useState(null);
  const [ showCreateFigure, setShowCreateFigure ] = useState(false);

  const classes = figureListStyles();

  useEffect( () => {
    json_get(
      `/api/${project_name}/figures`
    ).then(
      ({figures}) => setFigures(figures)
    )
  }, []);


  return (
    <main className={classes.figures}>
      <Grid container direction='column'>
        <Grid item container className={classes.title}><Typography variant='h5'>{project_name}</Typography></Grid>
        {
          figures ? figures.map(
            (figure,i) => <Figure figure={figure}/>
          ) : null
        }
        <Grid><Button variant='contained' onClick={() => setShowCreateFigure(true)} variant="text">Create Figure</Button></Grid>
        <CreateFigure
          project_name={project_name}
          open={showCreateFigure}
          handleClose={ () => setShowCreateFigure(false) }
        />
      </Grid>
    </main>
  );
}
