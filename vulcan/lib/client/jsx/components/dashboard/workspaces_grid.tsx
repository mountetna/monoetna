import React, {
  useCallback,
  useContext,
  useState,
  useEffect,
  useMemo
} from 'react';
import 'regenerator-runtime/runtime';

import {makeStyles} from '@material-ui/core/styles';

import {VulcanContext} from '../../contexts/vulcan_context';
import {Workflow, WorkspaceMinimal} from '../../api_types';
import WorkspaceCard from './cards/workspace_card';
import Grid from '@material-ui/core/Grid';
import { workflowByIdFromWorkflows } from '../../selectors/workflow_selectors';
import { updateWorkflowsWorkspaces } from '../../actions/vulcan_actions';

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
  allWorkflows,
  accessibleWorkspaces
}: {
  project_name: string;
  workflowId?: number | null;
  tags?: string[];
  searchString?: string;
  allWorkflows: Workflow[];
  accessibleWorkspaces: WorkspaceMinimal[];
}) {
  const {
    showErrors,
    updateWorkspace,
    deleteWorkspace,
    dispatch
  } = useContext(VulcanContext);

  if (!accessibleWorkspaces) return null;
  const [filteredWorkspaces, setFilteredWorkspaces] = useState<WorkspaceMinimal[]>(accessibleWorkspaces);

  const classes = useStyles();

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
    (workspace: WorkspaceMinimal) => {
      if (!workspace.workspace_id) return;
      const newTitle = prompt(
        'Please enter a new workspace title',
        workspace.name || ''
      );

      if (!newTitle) return;

      showErrors(
        updateWorkspace(project_name, workspace.workspace_id, newTitle, undefined)
        .then((updatedWorkspace) => {dispatch(updateWorkflowsWorkspaces())})
      );
      dispatch(updateWorkflowsWorkspaces());
    },
    [showErrors, updateWorkspace, project_name, accessibleWorkspaces]
  );

  const handleOnRemove = useCallback(
    (workspace: WorkspaceMinimal) => {
      if (!workspace.workspace_id) return;
      const doDelete = confirm(
        'This action will permanently delete this workspace and all of its files, a non-reversible action.'
      );
      if (!doDelete) return
      showErrors(deleteWorkspace(project_name, workspace.workspace_id))
      .then(() => {dispatch(updateWorkflowsWorkspaces())});
    },
    [dispatch, showErrors, deleteWorkspace, project_name, accessibleWorkspaces]
  );

  const handleUpdateTags = useCallback((tags: string[], workspace: WorkspaceMinimal) => {
    if (!workspace.workspace_id) return;
    // ToDo: Handle this update from the api return instead, to be sure we stay accurate
    showErrors(updateWorkspace(CONFIG.project_name, workspace.workspace_id as number, undefined, tags))
    .then(() => {dispatch(updateWorkflowsWorkspaces())});
  }, [CONFIG.project_name])

  const hasTag = useCallback(
    (workspace: WorkspaceMinimal) => {
      if (!tags || 0 === tags.length) return true;
      if (!workspace.tags || 0 === workspace.tags.length) return false;

      return (workspace.tags?.filter((t) => tags.includes(t)) || []).length > 0;
    },
    [tags]
  );

  const searchTerms = useMemo(() => {
    if (!searchString || '' === searchString) return '';
    // Break terms at spaces, unless "quoted together"
    let searchElements: string[] = [];
    let searchStringRemnant = searchString || ''
    const quotedRegex = /".+"/g
    if ((searchStringRemnant).match(quotedRegex)) {
      searchElements = (searchStringRemnant.match(quotedRegex) as string[]).map(s => s.replaceAll('"', ''))
      searchStringRemnant = searchStringRemnant.replaceAll(quotedRegex, '')
    }
    if (searchStringRemnant.trim().length>0) {
      searchElements.push(searchStringRemnant.trim().replaceAll('  ', ' ').replaceAll(' ', '|'))
    }
    return searchElements.join('|')
  }, [searchString])

  const matchesSearch = useCallback(
    (workspace: WorkspaceMinimal) => {
      if ('' === searchTerms) return true;

      const regex = new RegExp(searchTerms, 'i');
      const matchTo = [
        workspace.name,
        workspace.user_email,
        workspace.git_ref,
        workspace.git_sha,
        workflowByIdFromWorkflows(workspace.workflow_id, allWorkflows)?.name
      ]
      return matchTo.map(m => m?.match(regex)).some(m =>  m != null);
    },
    [searchTerms]
  );

  useEffect(() => {
    let results = accessibleWorkspaces
      .filter((workspace) => matchesSearch(workspace) && hasTag(workspace))
      .sort((a, b) => {
        if (null == a.workspace_id) return -1;
        if (null == b.workspace_id) return 1;
        return a.workspace_id - b.workspace_id;
      });
    // Prioritize published
    results = [
      ...results.filter((w) => (w.tags || []).includes('published')),
      ...results.filter((w) => !(w.tags || []).includes('published'))
    ]
    // Prioritize highlighted more
    setFilteredWorkspaces([
      ...results.filter((w) => (w.tags || []).includes('highlighted')),
      ...results.filter((w) => !(w.tags || []).includes('highlighted'))
    ]);
  }, [accessibleWorkspaces, matchesSearch, hasTag, workflowId]);

  return (
    <Grid className={classes.workspaces} container direction='row'>
      {filteredWorkspaces.map(
        (workspace: WorkspaceMinimal, index: number) => {
          return (
            <WorkspaceCard
              key={index}
              workspace={workspace}
              // onCopy={() => handleOnCopy(workspace)}
              onRemove={() => handleOnRemove(workspace)}
              onRename={() => handleOnRename(workspace)}
              updateTags={(tags: string[]) => handleUpdateTags(tags, workspace)}
            />
          );
        }
      )}
    </Grid>
  );
}
