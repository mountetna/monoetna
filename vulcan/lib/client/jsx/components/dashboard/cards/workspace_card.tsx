import React, {
  useCallback,
  useContext,
  useState,
  useEffect,
  useMemo
} from 'react';
import {makeStyles} from '@material-ui/core/styles';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {pushLocation} from 'etna-js/actions/location_actions';
import Menu from '@material-ui/core/Menu';
import MenuItem from '@material-ui/core/MenuItem';
import Card from '@material-ui/core/Card';
import CardHeader from '@material-ui/core/CardHeader';
import CardMedia from '@material-ui/core/CardMedia';
import CardContent from '@material-ui/core/CardContent';
import Avatar from '@material-ui/core/Avatar';
import Tooltip from '@material-ui/core/Tooltip';
import IconButton from '@material-ui/core/IconButton';
import MoreVertIcon from '@material-ui/icons/MoreVert';
import LaunchIcon from '@material-ui/icons/Launch';
import Grid from '@material-ui/core/Grid';

import {VulcanContext} from '../../../contexts/vulcan_context';
import {Workflow, Workspace, WorkspaceMinimal} from '../../../api_types';
import useUserHooks from '../../../contexts/useUserHooks';

import Tag from './tag';
import { workflowName } from '../../../selectors/workflow_selectors';
import Typography from '@material-ui/core/Typography';
import WorkspaceTagEditorModal from './tag_editor_modal';

const workspaceStyles = makeStyles((theme) => ({
  workspace: {
    border: '1px solid #eee',
    margin: '25px',
    width: '300px'
  },
  image: {
    cursor: 'pointer',
    borderTop: '1px solid #eee',
    borderBottom: '1px solid #eee'
  },
  defaultImage: {
    cursor: 'pointer',
    borderTop: '1px solid #eee',
    borderBottom: '1px solid #eee',
    opacity: 0.5,
    filter: 'grayscale(100%)'
  },
  author: {
    height: '30px',
    width: '30px',
    fontSize: '0.75rem'
  },
  title: {
    textOverflow: 'ellipsis',
    overflow: 'hidden',
    height: '2rem',
    whiteSpace: 'nowrap',
    maxWidth: '225px'
  },
  content: {
    width: '280px',
    padding: '10px 10px 0px 10px',
    '&:last-child': {
      paddingBottom: '14px'
    },
  },
  info: {
    maxWidth: '280px',
    padding: '0px'
  },
  tags: {
    maxWidth: '225px',
    maxHeight: '38px',
    overflow: 'hidden',
    position: 'relative',
    paddingTop: '5px',
    '&:after': {
      content: "''",
      position: 'absolute',
      zIndex: 1,
      top: 0,
      left: 0,
      backgroundImage:
        'linear-gradient(to bottom, rgba(255,255,255, 0), rgba(255,255,255,0) 90%, rgba(255,255,255, 1) 100%)',
      width: '100%',
    },
    '&:hover': {
      overflow: 'visible',
      zIndex: 2
    },
    '&:hover::after': {
      content: '',
      backgroundImage: 'none'
    }
  },
  tag: {
    marginBottom: '5px'
  }
}));

const authorInitials = ({user_email}: WorkspaceMinimal) => {
  if (!user_email) return '';

  let names = user_email.split("@")[0].split(/[^a-zA-Z]+/).map((n: string) => n.toUpperCase());
  return names.length > 1
    ? [names[0][0], names[names.length - 1][0]].join('')
    : names[0][0];
};

const workspaceImage = (
  workflow: Workflow | null,
  workspace: WorkspaceMinimal
): [Boolean, string] =>
  workspace.thumbnails && workspace.thumbnails.length > 0
    ? [false, workspace.thumbnails[0]]
    : [
        true,
        `/images/${workflow && workflow.image ? workflow.image : 'default.png'}`
      ];

const WorkspaceCard = ({
  workspace,
  // onCopy,
  onRename,
  onRemove,
  updateTags
}: {
  workspace: WorkspaceMinimal;
  // onCopy: () => void;
  onRename: () => void;
  onRemove: () => void;
  updateTags: (tags: string[]) => void;
}) => {
  const invoke = useActionInvoker();
  let {state} = useContext(VulcanContext);

  const {workflows, projectName} = state;

  const {canEdit, superuser} = useUserHooks();

  const [tagModalOpen, setTagModalOpen] = useState(false);

  const workflow = workflows && workspace.workflow_id
    ? workflows.find((w) => w.id == workspace.workflow_id) || null
    : null;

  const classes = workspaceStyles();

  const visitWorkspace = useCallback(() => {
    invoke(
      pushLocation(
        `/${projectName}/workspace/${workspace.workspace_id}`
      )
    );
  }, [invoke, projectName, workspace]);

  const [menuAnchor, setMenuAnchor] = useState(
    null as HTMLButtonElement | null
  );

  const handleClose = () => {
    setMenuAnchor(null);
  };

  // const handleOnCopy = useCallback(() => {
  //   handleClose();
  //   onCopy();
  // }, [onCopy]);

  const handleOnRename = useCallback(() => {
    handleClose();
    onRename();
  }, [onRename]);

  const handleOnRemove = useCallback(() => {
    handleClose();
    onRemove();
  }, [onRemove]);

  const openTagEditor = () => {
    handleClose();
    setTagModalOpen(true);
  }
  const closeTagEditor = () => {
    handleClose();
    setTagModalOpen(false);
  }

  const editor = useMemo(() => canEdit(workspace), [
    workspace,
    canEdit
  ]);

  const [defaultImage, image] = workspaceImage(workflow, workspace);

  const ContentLine = ({text, tooltip}: {text: string, tooltip?: string}) => {
    if (!tooltip) {
      tooltip = text;
    }
    return ( <Tooltip title={tooltip}>
      <Typography
        variant='body2'
        component="p"
        noWrap
        className={classes.info}
      >
        {text}
      </Typography>
    </Tooltip>)
  }

  return (
    <Card className={classes.workspace}>
      <Menu
        id={`${workspace.workspace_id}-dropdown-menu`}
        open={Boolean(menuAnchor)}
        anchorEl={menuAnchor}
        onClose={handleClose}
      >
        {/* <MenuItem onClick={handleOnCopy}>Copy</MenuItem> */}
        {editor ? <MenuItem onClick={handleOnRename}>Rename</MenuItem> : null}
        {editor ? <MenuItem onClick={handleOnRemove}>Remove</MenuItem> : null}
        {editor ? <MenuItem onClick={openTagEditor}>Edit Tags</MenuItem> : null}
      </Menu>
      <CardHeader
        title={workspace.name}
        titleTypographyProps={{
          variant: 'subtitle1',
          title: workspace.name,
          className: classes.title
        }}
        action={
          <div>
            <IconButton
              onClick={visitWorkspace}
              title={`Jump into this workspace!`}
              color={editor? 'primary' : 'default'}
            >
              <LaunchIcon/>
            </IconButton>
            {editor && <IconButton onClick={(e) => setMenuAnchor(e.currentTarget)}>
              <MoreVertIcon />
            </IconButton>}
          </div>
        }
      />
      <CardMedia
        className={defaultImage ? classes.defaultImage : classes.image}
        onClick={visitWorkspace}
        component='img'
        height='200'
        image={image}
        title={`Jump into this workspace!`}
      />
      <CardContent className={classes.content}>
        <Grid
          alignItems='flex-start'
          justifyContent='space-between'
          container
          direction='column'
          spacing={1}
        >
          <Grid item>
            <ContentLine text={`Pipeline: ${workflowName(workflow) || ''}`}/>
            <ContentLine
              text={
                workspace.git_ref != workspace.git_sha ?
                  `Version: ${workspace.git_ref} (${workspace.git_sha})` :
                  `Version: ${workspace.git_ref}`
              }
              tooltip={
                `Version Requested: ${workspace.git_ref} (Retrieved SHA: ${workspace.git_sha})`
              }
            />
            <ContentLine
              text={`Created: ${workspace.created_at.split(' ')[0]}`}
              tooltip={`Created: ${workspace.created_at.split(' ')[0]} (Id: ${workspace.workspace_id})`}
            />
          </Grid>
          <Grid
            item
            alignItems='center'
            justifyContent='space-between'
            container
            direction='row'
          >
            <Grid item xs={2}>
              <Tooltip title={workspace.user_email || ''}>
                <Avatar className={classes.author}>
                  {authorInitials(workspace)}
                </Avatar>
              </Tooltip>
            </Grid>
            <Grid item xs={10} container className={classes.tags}>
              {(workspace.tags || []).map((tag, index) => (
                <Tag className={classes.tag} key={index} label={tag} />
              ))}
            </Grid>
          </Grid>
        </Grid>
      </CardContent>
      <WorkspaceTagEditorModal workspace={workspace} open={tagModalOpen} onClose={closeTagEditor} handleUpdate={updateTags} showPushlishedOption/>
    </Card>
  );
};

export default WorkspaceCard;
