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

import {VulcanContext} from '../../contexts/vulcan_context';
import {VulcanFigureSession} from '../../api_types';
import useUserHooks from '../useUserHooks';

const figureStyles = makeStyles((theme) => ({
  figure: {
    border: '1px solid #eee'
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
    width: '250px',
    height: '2rem'
  }
}));

const authorInitials = ({author}: VulcanFigureSession) => {
  if (!author) return '';

  let names = author.split(/\s+/).filter((n: string) => /^[A-Z]/.test(n));
  return names.length > 1
    ? [names[0][0], names[names.length - 1][0]].join('')
    : names[0][0];
};

const figureImage = (workflow, figure) => figure.thumbnails.length > 0 ? figure.thumbnails[0] : `/images/${workflow ? workflow.image : 'default.png'}`;

const Figure = ({
  figureSession,
  onCopy,
  onRename,
  onRemove
}: {
  figureSession: VulcanFigureSession;
  onCopy: () => void;
  onRename: () => void;
  onRemove: () => void;
}) => {
  const invoke = useActionInvoker();
  let {state} = useContext(VulcanContext);
  const {workflows} = state;

  const {canEdit} = useUserHooks();

  const workflow = workflows
    ? workflows.find((w) => w.name == figureSession.workflow_name)
    : null;

  const classes = figureStyles();

  const visitFigure = useCallback(() => {
    invoke(
      pushLocation(
        `/${figureSession.project_name}/figure/${figureSession.figure_id}`
      )
    );
  }, [invoke, figureSession]);

  const [menuAnchor, setMenuAnchor] = useState(
    null as HTMLButtonElement | null
  );

  const handleClose = () => {
    setMenuAnchor(null);
  };

  const handleOnCopy = useCallback(() => {
    handleClose();
    onCopy();
  }, [onCopy]);

  const handleOnRename = useCallback(() => {
    handleClose();
    onRename();
  }, [onRename]);

  const handleOnRemove = useCallback(() => {
    handleClose();
    onRemove();
  }, [onRemove]);

  const editor = useMemo(() => canEdit(figureSession), [
    figureSession,
    canEdit
  ]);

  const image = figureImage(workflow, figureSession);

  return (
    <Card className={classes.figure}>
      <Menu
        id={`${figureSession.figure_id}-dropdown-menu`}
        open={Boolean(menuAnchor)}
        anchorEl={menuAnchor}
        onClose={handleClose}
      >
        <MenuItem onClick={handleOnCopy}>Copy</MenuItem>
        {editor ? <MenuItem onClick={handleOnRename}>Rename</MenuItem> : null}
        {editor ? <MenuItem onClick={handleOnRemove}>Remove</MenuItem> : null}
      </Menu>
      <CardHeader
        title={figureSession.title}
        titleTypographyProps={{variant: 'h6', className: classes.title}}
        subheader={figureSession.workflow_name.replace('.cwl', '')}
        action={
          <IconButton onClick={(e) => setMenuAnchor(e.currentTarget)}>
            <MoreVertIcon />
          </IconButton>
        }
      />
      <CardMedia
        className={
          image == 'default.png' ? classes.defaultImage : classes.image
        }
        onClick={visitFigure}
        component='img'
        height='140'
        image={image}
        title={figureSession.title || ''}
      />
      <CardContent>
        <Tooltip title={figureSession.author || ''}>
          <Avatar className={classes.author}>
            {authorInitials(figureSession)}
          </Avatar>
        </Tooltip>
      </CardContent>
    </Card>
  );
};

export default Figure;
