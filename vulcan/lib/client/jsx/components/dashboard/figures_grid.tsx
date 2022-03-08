import React, {
  useCallback,
  useContext,
  useState,
  useEffect,
  useMemo
} from 'react';
import 'regenerator-runtime/runtime';

import ImageList from '@material-ui/core/ImageList';
import ImageListItem from '@material-ui/core/ImageListItem';

import {makeStyles} from '@material-ui/core/styles';

import {VulcanContext} from '../../contexts/vulcan_context';
import {VulcanFigureSession} from '../../api_types';
import FigureCard from './figure';

const useStyles = makeStyles((theme) => ({
  title: {
    padding: '10px 0px',
    color: '#444'
  },
  figures: {
    padding: '15px',
    width: '100vw'
  }
}));

export default function FiguresTable({
  project_name,
  workflowName
}: {
  project_name: string;
  workflowName?: string;
}) {
  const {
    showErrors,
    fetchFigures,
    createFigure,
    updateFigure,
    deleteFigure,
    state
  } = useContext(VulcanContext);

  const [allFigureSessions, setAllFigureSessions] = useState<
    VulcanFigureSession[]
  >([]);
  const [filteredFigureSessions, setFilteredFigureSessions] = useState<
    VulcanFigureSession[]
  >([]);

  const classes = useStyles();

  useEffect(() => {
    showErrors(
      fetchFigures(
        project_name
      ).then(({figures}: {figures: VulcanFigureSession[]}) =>
        setAllFigureSessions(figures)
      )
    );
  }, [showErrors, fetchFigures, project_name]);

  useEffect(() => {
    if (workflowName) {
      setFilteredFigureSessions(
        allFigureSessions.filter(
          (figure) => figure.workflow_name === workflowName
        )
      );
    } else {
      setFilteredFigureSessions([...allFigureSessions]);
    }
  }, [workflowName, allFigureSessions]);

  const handleOnCopy = useCallback(
    (figure: VulcanFigureSession) => {
      const copy = {
        ...figure,
        figure_id: null,
        title: `${figure.title} - copy`
      };
      showErrors(
        createFigure(project_name, copy).then((newFigure) => {
          setAllFigureSessions([...allFigureSessions].concat([newFigure]));
        })
      );
    },
    [showErrors, createFigure, project_name, allFigureSessions]
  );

  const handleOnRename = useCallback(
    (figure: VulcanFigureSession) => {
      const newTitle = prompt(
        'Please enter a new figure title',
        figure.title || ''
      );

      if (!newTitle) return;

      showErrors(
        updateFigure(project_name, {
          ...figure,
          title: newTitle
        }).then((updatedFigure) => {
          const updated = allFigureSessions.map((oldFigure) => {
            if (oldFigure.figure_id === updatedFigure.figure_id) {
              return updatedFigure;
            }

            return oldFigure;
          });

          setAllFigureSessions(updated);
        })
      );
    },
    [showErrors, updateFigure, project_name, allFigureSessions]
  );

  const handleOnRemove = useCallback(
    (figure: VulcanFigureSession) => {
      if (!figure.figure_id) return;
      showErrors(
        deleteFigure(project_name, figure.figure_id).then(() => {
          const updated = allFigureSessions.filter((oldFigure) => {
            return oldFigure.figure_id !== figure.figure_id;
          });
          setAllFigureSessions(updated);
        })
      );
    },
    [showErrors, deleteFigure, project_name, allFigureSessions]
  );

  const filteredFigures = useMemo(() => {
    return filteredFigureSessions.sort((a, b) => {
      if (null == a.figure_id) return -1;
      if (null == b.figure_id) return 1;
      return a.figure_id - b.figure_id;
    });
  }, [filteredFigureSessions]);

  return (
    <ImageList cols={5} gap={30} rowHeight={330} className={classes.figures}>
      {filteredFigures.map((figure: VulcanFigureSession, index: number) => {
        return (
          <ImageListItem key={index}>
            <FigureCard
              figureSession={figure}
              onCopy={() => handleOnCopy(figure)}
              onRemove={() => handleOnRemove(figure)}
              onRename={() => handleOnRename(figure)}
            />
          </ImageListItem>
        );
      })}
    </ImageList>
  );
}
