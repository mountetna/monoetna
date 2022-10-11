import React, { useState, useEffect, useCallback } from 'react';
import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
import {makeStyles} from '@material-ui/core/styles';

const useStyles = makeStyles( theme => ({
  header: {
    padding: '10px 15px 5px'
  }
}));
const ProjectHeader = ({project_name, project_name_full, className, children }) => {
  const classes = useStyles();
  return <Grid item container className={`${classes.header} ${className}`} >
    <Typography color='primary' variant='h5'>
      {project_name_full || project_name}
    </Typography>
    {
      children
    }
  </Grid>
}

export default ProjectHeader;
