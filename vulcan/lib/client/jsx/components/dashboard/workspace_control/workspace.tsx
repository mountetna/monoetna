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
import Grid from '@material-ui/core/Grid';

import {VulcanContext} from '../../../contexts/vulcan_context';
import {Workflow, Workspace} from '../../../api_types';
import useUserHooks from '../../useUserHooks';

import Tag from '../tag';
import { workflowName } from '../../../selectors/workflow_selectors';

const workspaceStyles = makeStyles((theme) => ({
  content: {
    width: '100%',
    height: '100%'
  },
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
  tags: {
    maxWidth: '225px',
    maxHeight: '65px',
    overflow: 'hidden',
    position: 'relative',
    '&:after': {
      content: "''",
      position: 'absolute',
      zIndex: 1,
      top: 0,
      left: 0,
      backgroundImage:
        'linear-gradient(to bottom, rgba(255,255,255, 0), rgba(255,255,255,0) 90%, rgba(255,255,255, 1) 100%)',
      width: '100%',
      height: '65px'
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

const authorInitials = ({author}: Workspace) => {
  if (!author) return '';

  let names = author.split(/\s+/).filter((n: string) => /^[A-Z]/.test(n));
  return names.length > 1
    ? [names[0][0], names[names.length - 1][0]].join('')
    : names[0][0];
};

const workspaceImage = (
  workflow: Workflow | null,
  workspace: Workspace
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
  onRemove
}: {
  workspace: Workspace;
  // onCopy: () => void;
  onRename: () => void;
  onRemove: () => void;
}) => {
  const invoke = useActionInvoker();
  let {state} = useContext(VulcanContext);
  const {workflows, projectName} = state;

  const {canEdit} = useUserHooks();

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

  const editor = useMemo(() => canEdit(workspace), [
    workspace,
    canEdit
  ]);

  const [defaultImage, image] = workspaceImage(workflow, workspace);

  const title = workspace.name ? workspace.name : workflow ? `unnamed-${workflow.name}` : 'unnamed-workspace'

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
      </Menu>
      <CardHeader
        title={title}
        titleTypographyProps={{
          variant: 'subtitle1',
          title: title,
          className: classes.title
        }}
        subheader={workflowName(workflow) || ''}
        subheaderTypographyProps={{variant: 'subtitle2'}}
        action={
          <IconButton onClick={(e) => setMenuAnchor(e.currentTarget)}>
            <MoreVertIcon />
          </IconButton>
        }
      />
      <CardMedia
        className={defaultImage ? classes.defaultImage : classes.image}
        onClick={visitWorkspace}
        component='img'
        height='200'
        image={image}
        title={title}
      />
      <CardContent style={{width: '280px', height: '65px', padding: '10px'}}>
        <Grid
          alignItems='center'
          justifyContent='space-between'
          container
          className={classes.content}
        >
          <Grid item xs={1}>
            <Tooltip title={workspace.author || ''}>
              <Avatar className={classes.author}>
                {authorInitials(workspace)}
              </Avatar>
            </Tooltip>
          </Grid>
          <Grid item xs={11} container className={classes.tags}>
            {(workspace.tags || []).map((tag, index) => (
              <Tag className={classes.tag} key={index} label={tag} />
            ))}
          </Grid>
        </Grid>
      </CardContent>
    </Card>
  );
};

export default WorkspaceCard;
