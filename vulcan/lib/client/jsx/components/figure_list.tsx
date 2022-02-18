import React, {
  useCallback,
  useContext,
  useState,
  useEffect,
  useMemo
} from 'react';
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
import {VulcanFigure, VulcanFigureSession, Workflow} from '../api_types';
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

  const projectWorkflows = useMemo(() => {
    return workflows
      ? workflows.filter(
          ({projects}) => projects && projects.includes(project_name)
        )
      : [];
  }, [workflows, project_name]);

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

  const numColumns = useMemo(() => {
    return Math.min(projectWorkflows.length, 4);
  }, [projectWorkflows]);

  return (
    <Dialog onClose={close} open={open} maxWidth='xl'>
      <DialogTitle>Select a workflow</DialogTitle>
      <DialogContent dividers>
        {projectWorkflows.length ? (
          <ImageList rowHeight='auto' cols={numColumns}>
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
  const {
    showErrors,
    fetchFigures,
    createFigure,
    updateFigure,
    deleteFigure
  } = useContext(VulcanContext);

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

  const handleOnCopy = useCallback(
    (figure: VulcanFigureSession) => {
      const copy = {
        ...figure,
        figure_id: null,
        title: `${figure.title} - copy`
      };
      showErrors(
        createFigure(project_name, copy).then((newFigure) => {
          setFigureSessions([...figureSessions].concat([newFigure]));
        })
      );
    },
    [showErrors, createFigure, project_name, figureSessions]
  );

  const handleOnRename = useCallback(
    (figure: VulcanFigureSession) => {
      const newTitle = prompt(
        'Please enter a new figure title',
        figure.title || ''
      );
      showErrors(
        updateFigure(project_name, {
          ...figure,
          title: newTitle
        }).then((updatedFigure) => {
          const updated = figureSessions.map((oldFigure) => {
            if (oldFigure.figure_id === updatedFigure.figure_id) {
              return updatedFigure;
            }

            return oldFigure;
          });

          setFigureSessions(updated);
        })
      );
    },
    [showErrors, updateFigure, project_name, figureSessions]
  );

  const handleOnRemove = useCallback(
    (figure: VulcanFigureSession) => {
      if (!figure.figure_id) return;

      showErrors(
        deleteFigure(project_name, figure.figure_id).then(() => {
          const updated = figureSessions.filter((oldFigure) => {
            return oldFigure.figure_id !== figure.figure_id;
          });
          setFigureSessions(updated);
        })
      );
    },
    [showErrors, deleteFigure, project_name, figureSessions]
  );

  return (
    <main className={classes.figures}>
      <Grid container direction='column'>
        <Grid item className={classes.title}>
          <Typography variant='h5'>{project_name}</Typography>
        </Grid>
        <Grid container direction='row'>
          {figureSessions
            ? figureSessions.map((figureSession, i) => (
                <Figure
                  key={i}
                  figureSession={figureSession}
                  onCopy={() => handleOnCopy(figureSession)}
                  onRename={() => handleOnRename(figureSession)}
                  onRemove={() => handleOnRemove(figureSession)}
                />
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
