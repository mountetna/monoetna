import React, { useState, useRef } from 'react';
require('./Nav.css');
import { isSuperuser, isSuperEditor, isSuperViewer } from '../utils/janus';
import AppBar from '@material-ui/core/AppBar';
import Toolbar from '@material-ui/core/Toolbar';
import IconButton from '@material-ui/core/IconButton';
import AppsIcon from '@material-ui/icons/Apps';
import ClickAwayListener from '@material-ui/core/ClickAwayListener';
import Grow from '@material-ui/core/Grow';
import Paper from '@material-ui/core/Paper';
import Popper from '@material-ui/core/Popper';
import Grid from '@material-ui/core/Grid';
import Icon from './icon';

const ICONS = {
  superuser: 'user-ninja',
  supereditor: 'user-nurse',
  superviewer: 'user-tie',
  administrator: 'user-astronaut',
  editor: 'user-edit',
  viewer: 'user'
};

const APPS = [
  { name: "timur", description: "Model" },
  { name: "metis", description: "Store" },
  { name: "janus", description: "Auth" },
  { name: "polyphemus", description: "Link" },
  { name: "vulcan", description: "Analyze" },
  { name: "gnomon", description: "Name" },
]

const Login = ({ user }) => {
  if (!user) return null;

  let { name, permissions } = user;

  let role = (permissions[CONFIG.project_name] || {}).role;

  if (isSuperuser(user)) role = 'superuser';
  else if (isSuperEditor(user)) role = 'supereditor';
  else if (isSuperViewer(user)) role = 'superviewer';

  return (
    <div className='etna-login' title={role}>
      <Icon className='etna-user' icon={ICONS[role]} />
      {name}
    </div>
  );
};

const Logo = ({ LogoImage, app }) =>
  <div className='etna-logo'>
    <a href='/' title={titelize(app)}>
      <LogoImage />
    </a>
  </div>;

const Link = ({ app }) => {
  const name = app.name
  const host_key = `${name}_host`;

  const link = CONFIG[host_key] ? new URL(...[CONFIG.project_name, CONFIG[host_key]].filter(_ => _)) : undefined;

  return <a href={link} className='etna-link' title={titelize(name)}>
    <img src={`/images/${name}.svg`} />
    <div className="description">{app.description}</div>
  </a>;
}

const AppsMenu = ({ currentApp }) => {
  const [open, setOpen] = useState(false);
  const anchorRef = useRef(null);

  const handleToggle = () => {
    setOpen((prevOpen) => !prevOpen);
  };

  const handleClose = (event) => {
    setOpen(false);
  };

  return (
    <div className='apps-menu-container'>
      <IconButton
        ref={anchorRef}
        onClick={handleToggle}
        aria-label="Show Apps"
        aria-haspopup="true"
        aria-controls={open ? 'apps-menu-list' : undefined}
        disableRipple
        disableFocusRipple
      >
        <AppsIcon />
      </IconButton>
      <Popper
        open={open}
        anchorEl={anchorRef.current}
        placement='bottom-end'
        role={undefined}
        transition
        disablePortal
      >
        {({ TransitionProps }) => (
          <Grow
            {...TransitionProps}
            style={{ transformOrigin: 'right top' }}
          >
            <Paper className='apps-menu' variant='outlined'>
              <ClickAwayListener onClickAway={handleClose}>
                <Grid container id="apps-menu-list">
                  {
                    APPS.map(app => {
                      return (
                        <Grid item xs={6} key={app.name} className='apps-menu-item'>
                          <Link app={app} />
                        </Grid>
                      )
                    })
                  }
                </Grid>
              </ClickAwayListener>
            </Paper>
          </Grow>
        )}
      </Popper>
    </div>
  );
}

function findValidChildren(children) {
  if (children) {
    // Only return non-null children from an Array, or
    //  the singular child if provided. When only
    //  one child is passed in, it is an Object, not an Array,
    //  and has no `filter` method.
    if (Array.isArray(children)) {
      if (children.filter(_ => _).length) {
        return children;
      }
    } else {
      return children;
    }
  }
  return <div style={{ flex: 1 }} />;
}

function titelize(word) {
  return word.charAt(0).toUpperCase() + word.slice(1)
}

// TODO
// make responsive?
const Nav = ({ logo, app, children, user }) => {
  return (
    <AppBar position="sticky" className="etna-nav">
      <Toolbar disableGutters className="etna-nav-toolbar">
        <Logo LogoImage={logo} app={app} />
        {findValidChildren(children)}
        <Login user={user} />
        <AppsMenu currentApp={app} />
      </Toolbar>
    </AppBar>
  );
}



export default Nav;
