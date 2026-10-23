import React, {
  useCallback,
  useEffect,
  useState,
  useMemo,
  useContext
} from 'react';

import Dialog from '@material-ui/core/Dialog';
import DialogTitle from '@material-ui/core/DialogTitle';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import Autocomplete from '@material-ui/lab/Autocomplete';
import Button from '@material-ui/core/Button';
import {makeStyles} from '@material-ui/core/styles';
import TextField from '@material-ui/core/TextField';
import Tag from './tag';
import { Workspace, WorkspaceMinimal } from '../../../api_types';
import * as _ from 'lodash';
import {VulcanContext} from '../../../contexts/vulcan_context';
import useUserHooks from '../../../contexts/useUserHooks';


const useStyles = makeStyles((theme) => ({
  title: {
    '&:last-child': {flex: '1 1 auto'}
  },
  titleText: {
    textOverflow: 'ellipsis',
    overflow: 'hidden'
  },
  tags: {
    padding: '12.5px !important'
  },
  editTags: {
    width: '600px'
  }
}));

export default function WorkspaceTagEditorModal({
  workspace,
  open,
  onClose,
  handleUpdate,
  tagOptions,
  showPushlishedOption = false,
}: {
  workspace: Workspace | WorkspaceMinimal;
  open: boolean;
  onClose: () => void;
  handleUpdate: (tags: string[]) => void;
  tagOptions?: string[]
  showPushlishedOption?: boolean
}) {

  const classes = useStyles();

  let {state} = useContext(VulcanContext);
  const {workflows, workspaces} = useMemo(() => {
    return state
  }, [state.workflows, state.workspaces])
  const {canEdit} = useUserHooks();
  const userVisibleWorkspaceTags = useMemo(() => {
    const projectWorkflowIds = workflows.filter(w => w.project_name == CONFIG.project_name).map(w => w.id)
    const projectWorkspaces = workspaces.filter((w) => projectWorkflowIds.includes(w.workflow_id))
    return [...new Set(projectWorkspaces.filter((w) => canEdit(w) || (w.tags || []).includes('published')).map(w => w.tags || []).flat())]
  }, [workflows, workspaces, CONFIG.project_name]);

  const [localTags, setLocalTags] = useState<string[]>(workspace.tags || []);

  function resetTags() {
    setLocalTags(workspace.tags || []);
  }

  useEffect(() => {
    resetTags();
  }, [workspace.tags, open])

  const tapOptionsUse = useMemo(() => {
    const tagOptionsUse = tagOptions || userVisibleWorkspaceTags || []
    if (!tagOptionsUse?.includes('published') && showPushlishedOption) {
      tagOptionsUse?.push('published')
    }
    if (!tagOptionsUse?.includes('highlighted') && showPushlishedOption) {
      tagOptionsUse?.push('highlighted')
    }
    return tagOptionsUse
  }, [tagOptions, userVisibleWorkspaceTags])

  // const isPublic = useMemo(() => (workspace.tags || []).includes('published'), [workspace.tags]);
  // function togglePublicTag() {
  //   const currentTags = workspace.tags || [];
  //   const newTags = currentTags.includes('published') ?
  //   currentTags.filter((v)=>v!='published') :
  //     [...currentTags, 'published']
  //   setLocalTags(newTags);
  //   handleUpdate(newTags);
  // }

  return (
    <>
      <Dialog
        maxWidth='md'
        open={open}
        onClose={onClose}
      >
        <DialogTitle id='tag-editor'>Edit Tags</DialogTitle>
        <DialogContent className={classes.editTags}>
        <Autocomplete
          fullWidth
          multiple
          freeSolo
          className='figure-edit-tag-autocomplete'
          classes={{
            input: classes.tags
          }}
          value={localTags}
          id='figure-edit-tags-filter'
          options={tapOptionsUse}
          renderInput={(params: any) => (
            <TextField {...params} label='Tags' variant='outlined' />
          )}
          renderTags={(tags: string[], getTagProps: any) =>
            tags.map((tag, index) => (
              <Tag {...getTagProps({index})} label={tag} />
            ))
          }
          renderOption={(option: string, state: any) => (
            <span>{option}</span>
          )}
          filterOptions={(options: string[], state: any) => {
            let regex = new RegExp(state.inputValue);
            return options.filter((o) => regex.test(o));
          }}
          onChange={(e: any, v: string[]) => setLocalTags(v)}
        />
        </DialogContent>
        <DialogActions>
        <Button
          onClick={() => {
            handleUpdate(localTags)
            onClose()
          }}
          color='primary'
          disabled={_.isEqual(localTags,workspace.tags)}
        >
          Save Tags
        </Button>
        <Button onClick={onClose} color='primary'>
          Close
        </Button>
        </DialogActions>
      </Dialog>
    </>
  )
}
