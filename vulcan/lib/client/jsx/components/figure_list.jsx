import React, {useCallback, useContext, useState, useEffect} from 'react';
import 'regenerator-runtime/runtime';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {pushLocation} from 'etna-js/actions/location_actions';

import {VulcanContext} from '../contexts/vulcan_context';
import WorkflowCard from '../components/dashboard/card';
import Figure from '../components/figure';

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
import {json_get} from 'etna-js/utils/fetch';

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

const CreateFigure = ({handleClose, open, project_name}) => {
  let {state} = useContext(VulcanContext);
  const {workflows} = state;

  const [selectedWorkflow, selectWorkflow] = useState(null);
  const classes = createFigureStyles();

  const projectWorkflows = workflows
    ? workflows.filter(({projects}) => projects.includes(project_name))
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
    [invoke]
  );

  const close = useCallback(() => {
    handleClose();
    selectWorkflow(null);
  }, [handleClose]);

  return (
    <Dialog onClose={close} open={open} maxWidth='xl'>
      <DialogTitle>Select a workflow</DialogTitle>
      <DialogContent dividers>
        {projectWorkflows.length ? (
          <ImageList rowHeight='auto' cols={3}>
            {projectWorkflows.map((w, ind) => (
              <ImageListItem key={ind}>
                <WorkflowCard
                  className={selectedWorkflow == w ? classes.selected : null}
                  workflow={w}
                  key={ind}
                  onClick={() =>
                    selectWorkflow(selectedWorkflow == w ? null : w)
                  }
                />
              </ImageListItem>
            ))}
          </ImageList>
        ) : (
          // <Grid className={classes.workflows} container direction='column'>
          //   {projectWorkflows.map((w, ind) => (
          //     <WorkflowCard
          //       className={selectedWorkflow == w ? classes.selected : null}
          //       workflow={w}
          //       key={ind}
          //       onClick={() => selectWorkflow(selectedWorkflow == w ? null : w)}
          //     />
          //   ))}
          // </Grid>
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

export default function FigureList({project_name}) {
  const invoke = useActionInvoker();
  const [figures, setFigures] = useState(null);
  const [showCreateFigure, setShowCreateFigure] = useState(false);

  const classes = figureListStyles();

  useEffect(() => {
    json_get(`/api/${project_name}/figures`).then(({figures}) =>
      setFigures(figures)
    );
  }, []);

  return (
    <main className={classes.figures}>
      <Grid container direction='column'>
        <Grid item className={classes.title}>
          <Typography variant='h5'>{project_name}</Typography>
        </Grid>
        <Grid container direction='row'>
          {figures
            ? figures.map((figure, i) => <Figure key={i} figure={figure} />)
            : null}
        </Grid>
        <Grid>
          <Button
            variant='contained'
            onClick={() => setShowCreateFigure(true)}
            variant='text'
          >
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
