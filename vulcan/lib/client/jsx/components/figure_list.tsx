import React, {useCallback, useContext, useState, useEffect} from 'react';
import 'regenerator-runtime/runtime';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {pushLocation} from 'etna-js/actions/location_actions';

import {VulcanContext} from '../contexts/vulcan_context';
import WorkflowCard from './dashboard/card';
import Figure from './figure';

import ImageList from '@material-ui/core/ImageList';
import ImageListItem from '@material-ui/core/ImageListItem';
import Grid from '@material-ui/core/Grid';
import Dialog from '@material-ui/core/Dialog';
import DialogContent from '@material-ui/core/DialogContent';
import DialogActions from '@material-ui/core/DialogActions';
import DialogTitle from '@material-ui/core/DialogTitle';
import Button from '@material-ui/core/Button';
import Typography from '@material-ui/core/Typography';

import {makeStyles} from '@material-ui/core/styles';
import {VulcanFigureSession, Workflow} from '../api_types';
import {showMessages} from 'etna-js/actions/message_actions';

const createFigureStyles = makeStyles((theme) => ({
  none: {
    color: '#f44'
  },
  workflows: {},
  selected: {
    border: '1px solid red',
    background: 'rgba(100,0,0,0.1)',
    boxShadow: '0 0 0 15px #fee, 0 0 4px 15px #faa'
  }
}));

const CreateFigure = ({
  handleClose,
  open,
  project_name
}: {
  open: boolean;
  project_name: string;
  handleClose: () => void;
}) => {
  let {state} = useContext(VulcanContext);
  const {workflows} = state;

  const [selectedWorkflow, setSelectedWorkflow] = useState(
    null as Workflow | null
  );
  const classes = createFigureStyles();

  const projectWorkflows = workflows
    ? workflows.filter(
        ({projects}) => projects && projects.includes(project_name)
      )
    : [];

  const invoke = useActionInvoker();

  const visitNewFigure = useCallback(
    (workflow) => {
      invoke(
        pushLocation(
          `/${project_name}/figure/new/${workflow.name.replace('.cwl', '')}`
        )
      );
    },
    [invoke, project_name]
  );

  const close = useCallback(() => {
    handleClose();
    setSelectedWorkflow(null);
  }, [handleClose]);

  return (
    <Dialog onClose={close} open={open} maxWidth='xl'>
      <DialogTitle>Select a workflow</DialogTitle>
      <DialogContent dividers>
        {projectWorkflows.length ? (
          <ImageList
            rowHeight='auto'
            cols={Math.min(projectWorkflows.length, 3)}
          >
            {projectWorkflows.map((w, ind) => (
              <ImageListItem key={ind}>
                <WorkflowCard
                  className={selectedWorkflow == w ? classes.selected : null}
                  workflow={w}
                  key={ind}
                  onClick={() =>
                    setSelectedWorkflow(selectedWorkflow == w ? null : w)
                  }
                />
              </ImageListItem>
            ))}
          </ImageList>
        ) : (
          <Grid item className={classes.none}>
            <em>No workflows</em>
          </Grid>
        )}
      </DialogContent>
      <DialogActions>
        {selectedWorkflow && (
          <Typography variant='body1'>
            {selectedWorkflow.displayName}
          </Typography>
        )}
        <Button
          variant='contained'
          disabled={!selectedWorkflow}
          onClick={() => visitNewFigure(selectedWorkflow)}
        >
          Create Figure
        </Button>
      </DialogActions>
    </Dialog>
  );
};

const figureListStyles = makeStyles((theme) => ({
  title: {
    padding: '10px 0px',
    color: '#444'
  },
  figures: {
    padding: '15px'
  }
}));

export default function FigureList({project_name}: {project_name: string}) {
  const {showErrors, fetchFigures} = useContext(VulcanContext);

  const [figureSessions, setFigureSessions] = useState(
    [] as VulcanFigureSession[]
  );
  const [showCreateFigure, setShowCreateFigure] = useState(false);

  const classes = figureListStyles();

  useEffect(() => {
    showErrors(
      fetchFigures(project_name).then(({figures}) => setFigureSessions(figures))
    );
  }, [showErrors, fetchFigures, project_name]);

  return (
    <main className={classes.figures}>
      <Grid container direction='column'>
        <Grid item className={classes.title}>
          <Typography variant='h5'>{project_name}</Typography>
        </Grid>
        <Grid container direction='row'>
          {figureSessions
            ? figureSessions.map((figureSession, i) => (
                <Figure key={i} figureSession={figureSession} />
              ))
            : null}
        </Grid>
        <Grid>
          <Button onClick={() => setShowCreateFigure(true)} variant='text'>
            Create Figure
          </Button>
        </Grid>
        <CreateFigure
          project_name={project_name}
          open={showCreateFigure}
          handleClose={() => setShowCreateFigure(false)}
        />
      </Grid>
    </main>
  );
}
