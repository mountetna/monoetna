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
import TextField from '@material-ui/core/TextField';
import {makeStyles} from '@material-ui/core/styles';
import InputAdornment from '@material-ui/core/InputAdornment';
import SearchIcon from '@material-ui/icons/Search';
import Autocomplete from '@material-ui/lab/Autocomplete';

import {VulcanContext} from '../../contexts/vulcan_context';
import {VulcanFigureSession} from '../../api_types';
import FigureCard from './figure';
import Grid from '@material-ui/core/Grid';

const useStyles = makeStyles((theme) => ({
  figures: {
    boxShadow: '0 0 15px 0px #f5f5f5 inset'
  }
}));

export default function FiguresGrid({
  project_name,
  workflowName,
  tags,
  searchString,
  setSearchString,
  setTags
}: {
  project_name: string;
  workflowName?: string;
  tags?: string[];
  searchString?: string;
  setSearchString: Function;
  setTags: Function;
}) {
  const {
    showErrors,
    fetchFigures,
    createFigure,
    updateFigure,
    deleteFigure
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
  }, []);

  const handleOnCopy = useCallback(
    (figure: VulcanFigureSession) => {
      const copy = {
        ...figure,
        figure_id: null,
        title: `${figure.title} - copy`,
        tags: []
      };
      showErrors(
        createFigure(project_name, copy).then((newFigure) => {
          setAllFigureSessions([...allFigureSessions].concat([newFigure]));
          setTags([]);
          setSearchString(figure.title);
        })
      );
    },
    [
      showErrors,
      createFigure,
      project_name,
      allFigureSessions,
      setTags,
      setSearchString
    ]
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

  const hasTag = useCallback(
    (figure: VulcanFigureSession) => {
      if (!tags || 0 === tags.length) return true;
      if (!figure.tags || 0 === figure.tags.length) return false;

      return (figure.tags?.filter((t) => tags.includes(t)) || []).length > 0;
    },
    [tags]
  );

  const matchesSearch = useCallback(
    (figure: VulcanFigureSession) => {
      if (!searchString || '' === searchString) return true;

      const regex = new RegExp(searchString, 'i');
      return (
        figure.title?.match(regex) ||
        figure.author?.match(regex) ||
        figure.workflow_name?.match(regex)
      );
    },
    [searchString]
  );

  useEffect(() => {
    let results = allFigureSessions
      .filter((figure) => matchesSearch(figure) && hasTag(figure))
      .sort((a, b) => {
        if (null == a.figure_id) return -1;
        if (null == b.figure_id) return 1;
        return a.figure_id - b.figure_id;
      });

    if (workflowName) {
      setFilteredFigureSessions(
        results.filter((figure) => figure.workflow_name === workflowName)
      );
    } else {
      setFilteredFigureSessions([...results]);
    }
  }, [allFigureSessions, matchesSearch, hasTag, workflowName]);

  return (
    <Grid className={classes.figures} container direction='row'>
      {filteredFigureSessions.map(
        (figure: VulcanFigureSession, index: number) => {
          return (
            <FigureCard
              key={index}
              figureSession={figure}
              onCopy={() => handleOnCopy(figure)}
              onRemove={() => handleOnRemove(figure)}
              onRename={() => handleOnRename(figure)}
            />
          );
        }
      )}
    </Grid>
  );
}
