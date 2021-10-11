import React, { useState, useCallback } from 'react';
import Collapse from '@material-ui/core/Collapse';
import CardContent from '@material-ui/core/CardContent';
import CardActions from '@material-ui/core/CardActions';
import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
import {makeStyles} from '@material-ui/core/styles';

const useStyles = makeStyles( theme => ({
  pane: {
    padding: '8px',
  },
  paneheader: {
    height: '20px',
    alignItems: 'center'
  },
  header: {
    marginBottom: '5px',
    height: '35px'
  },
  title: {
    flexBasis: '200px'
  }
}));

export const EtlPaneHeader = ({title,children}:{title:string,children?:React.ReactNode}) => {
  const classes = useStyles();

  return <Grid className={classes.header} container alignItems='center'>
    <Grid item className={classes.title}>
      <Typography>{title}</Typography>
    </Grid>
    { children }
  </Grid>
}

const EtlPane = ({mode, selected, children }:{mode:string, selected:string|null,children:React.ReactNode}) => {
  const classes = useStyles();

  return <Collapse in={mode == selected} timeout='auto' unmountOnExit>
    <CardContent className={classes.pane}>
      { children }
    </CardContent>
  </Collapse>
};

export default EtlPane;
