import React, {useCallback, useContext, useState, useEffect} from 'react';
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


import {VulcanContext} from '../contexts/vulcan_context';

const figureStyles = makeStyles( theme => ({
  figure: {
    width: 350,
    marginRight: 20
  },
  image: {
    cursor: 'pointer',
    borderTop: '1px solid #eee',
    borderBottom: '1px solid #eee'
  },
  author: {
    height: '30px',
    width: '30px',
    fontSize: '0.75rem'
  }
}));

const authorInitials = ({author}) => {
  let names = author.split(/\s+/).filter(n => (/^[A-Z]/).test(n))
  return (names.length > 1) ? [ names[0][0], names[names.length-1][0] ].join('') : names[0][0];
}

const Figure = ({figure}) => {
  const invoke = useActionInvoker();
  let {state} = useContext(VulcanContext);
  const {workflows} = state;

  const workflow = workflows ? workflows.find( w => w.name == figure.workflow_name ) : null

  const classes = figureStyles();

  const visitFigure = useCallback((figure) => {
    invoke(pushLocation(`/${figure.project_name}/figure/${figure.figure_id}`));
  }, [invoke]);

  const [ menuAnchor, setMenuAnchor ] = useState(null);

  const handleClose = () => { setMenuAnchor(null); }

  return <Card className={classes.figure}>
    <Menu
      id="simple-menu"
      open={Boolean(menuAnchor)}
      anchorEl={ menuAnchor }
      onClose={ handleClose }
    >
        <MenuItem onClick={handleClose}>Copy</MenuItem>
        <MenuItem onClick={handleClose}>Rename</MenuItem>
        <MenuItem onClick={handleClose}>Remove</MenuItem>
    </Menu>
    <CardHeader
      title={figure.title}
      titleTypographyProps={ { variant: 'h6' } }
      subheader={figure.workflow_name.replace('.cwl','')}
      action={
        <IconButton onClick={ (e) => setMenuAnchor(e.currentTarget) }>
          <MoreVertIcon />
        </IconButton>
      }/>
    <CardMedia
      className={classes.image}
      onClick={ () => visitFigure(figure) }
      component="img"
      height="140"
      image={`/images/${workflow ? workflow.image : 'default.png'}`}
      title={figure.title}
    />
    <CardContent>
      <Tooltip title={ figure.author }>
        <Avatar className={ classes.author } >{ authorInitials(figure) }</Avatar>
      </Tooltip>
    </CardContent>
  </Card>
}

export default Figure;
