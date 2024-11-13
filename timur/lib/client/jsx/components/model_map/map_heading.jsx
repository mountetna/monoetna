import React from 'react';

import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
import {makeStyles} from '@material-ui/core/styles';

const headingStyles = makeStyles(theme => ({
  heading: {
    borderBottom: '1px solid rgba(34, 139, 34, 0.1)',
    flex: '0 0 48px',
    alignItems: 'center'
  },
  name: {
    color: 'gray'
  },
  title: {
    color: 'forestgreen'
  },
  buttons: {
    width: 'auto',
    flex: 1,
    alignItems: 'center',
    justifyContent: 'right'
  }
}));

const MapHeading = ({className, name, title, children}) => {
  const classes = headingStyles();
  return <Grid className={ `${classes.heading} ${className}` } container>
      <Typography variant='h6' className={classes.name}>{name}</Typography>
      &nbsp;
      <Typography variant='h6' className={classes.title}>{title}</Typography>

      <Grid className={classes.buttons} container>
      { children }
      </Grid>
    </Grid>;
};

export default MapHeading;
