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

import {VulcanContext} from '../../../contexts/vulcan_context';
import {WorkflowsResponse, Workspace, Workspaces} from '../../../api_types';
import WorkspaceCard from './workspace';
import Grid from '@material-ui/core/Grid';
import { workflowByIdFromWorkflows } from '../../../selectors/workflow_selectors';

const useStyles = makeStyles((theme) => ({
  workspaces: {
    boxShadow: '0 0 15px 0px #f5f5f5 inset'
  }
}));

export default function WorkspacesGrid({
  project_name,
  workflowId,
  tags,
  searchString,
  setSearchString,
  setTags
}: {
  project_name: string;
  workflowId?: number | null;
  tags?: string[];
  searchString?: string;
  setSearchString: Function;
  setTags: Function;
}) {
  const {
    showErrors,
    getWorkspaces,
    getWorkflows,
    createWorkspace,
    updateWorkspace,
    deleteWorkspace
  } = useContext(VulcanContext);

  const [allWorkflows, setAllWorkflows] = useState<WorkflowsResponse>([]);
  const [allWorkspaces, setAllWorkspaces] = useState<Workspaces>([]);
  const [filteredWorkspaces, setFilteredWorkspaces] = useState<Workspaces>([]);

  const classes = useStyles();

  useEffect(() => {
    showErrors(
      getWorkspaces(
        project_name
      ).then(({workspaces}: {workspaces: Workspaces}) =>
        setAllWorkspaces(workspaces)
      )
    );
    showErrors(
      getWorkflows(
        project_name
      ).then((workflows: WorkflowsResponse) =>
        setAllWorkflows(workflows)
      )
    );
  }, []);

  // const handleOnCopy = useCallback(
  //   (figure: VulcanFigureSession) => {
  //     const copy = {
  //       ...figure,
  //       figure_id: null,
  //       title: `${figure.title} - copy`,
  //       tags: []
  //     };
  //     showErrors(
  //       createFigure(project_name, copy).then((newFigure) => {
  //         setAllFigureSessions([...allFigureSessions].concat([newFigure]));
  //         setTags([]);
  //         setSearchString(figure.title);
  //       })
  //     );
  //   },
  //   [
  //     showErrors,
  //     createFigure,
  //     project_name,
  //     allFigureSessions,
  //     setTags,
  //     setSearchString
  //   ]
  // );

  const handleOnRename = useCallback(
    (workspace: Workspace) => {
      if (!workspace.workspace_id) return;
      const newTitle = prompt(
        'Please enter a new workspace title',
        workspace.name || ''
      );

      if (!newTitle) return;

      showErrors(
        updateWorkspace(project_name, workspace.workspace_id, newTitle, workspace.tags)
        .then((updatedWorkspace) => {
          const updated = allWorkspaces.map((oldWorkspace) => {
            if (oldWorkspace.workspace_id === updatedWorkspace.workspace_id) {
              return updatedWorkspace;
            }

            return oldWorkspace;
          });

          setAllWorkspaces(updated);
        })
      );
    },
    [showErrors, updateWorkspace, project_name, allWorkspaces]
  );

  const handleOnRemove = useCallback(
    (workspace: Workspace) => {
      if (!workspace.workspace_id) return;
      showErrors(
        deleteWorkspace(project_name, workspace.workspace_id).then(() => {
          const updated = allWorkspaces.filter((oldWorkspace) => {
            return oldWorkspace.workspace_id !== workspace.workspace_id;
          });
          setAllWorkspaces(updated);
        })
      );
    },
    [showErrors, deleteWorkspace, project_name, allWorkspaces]
  );

  const hasTag = useCallback(
    (workspace: Workspace) => {
      if (!tags || 0 === tags.length) return true;
      if (!workspace.tags || 0 === workspace.tags.length) return false;

      return (workspace.tags?.filter((t) => tags.includes(t)) || []).length > 0;
    },
    [tags]
  );

  const matchesSearch = useCallback(
    (workspace: Workspace) => {
      if (!searchString || '' === searchString) return true;

      const regex = new RegExp(searchString, 'i');
      return (
        workspace.name?.match(regex) ||
        workspace.author?.match(regex) ||
        workflowByIdFromWorkflows(workspace.workflow_id, allWorkflows)?.name.match(regex)
      );
    },
    [searchString]
  );

  useEffect(() => {
    let results = allWorkspaces
      .filter((workspace) => matchesSearch(workspace) && hasTag(workspace))
      .sort((a, b) => {
        if (null == a.workspace_id) return -1;
        if (null == b.workspace_id) return 1;
        return a.workspace_id - b.workspace_id;
      });

    if (workflowId) {
      setFilteredWorkspaces(
        results.filter((workspace) => workspace.workflow_id === workflowId)
      );
    } else {
      setFilteredWorkspaces([...results]);
    }
  }, [allWorkspaces, matchesSearch, hasTag, workflowId]);

  return (
    <Grid className={classes.workspaces} container direction='row'>
      {filteredWorkspaces.map(
        (workspace: Workspace, index: number) => {
          return (
            <WorkspaceCard
              key={index}
              workspace={workspace}
              // onCopy={() => handleOnCopy(workspace)}
              onRemove={() => handleOnRemove(workspace)}
              onRename={() => handleOnRename(workspace)}
            />
          );
        }
      )}
    </Grid>
  );
}
