import React, { useState, useEffect, useCallback } from 'react';
import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
import Link from '@material-ui/core/Link';
import {makeStyles} from '@material-ui/core/styles';

const useStyles = makeStyles( theme => ({
  header: {
    padding: '10px 15px 5px'
  },
  link: {
    '&:hover': {
      color: theme.palette.secondary.main
    }
  }
}));
const ProjectHeader = ({project_name, project_name_full, className, children, link=true }) => {
  const classes = useStyles();
  const title = <Typography variant='h5'>
    {project_name_full || project_name}
  </Typography>;
  const projectPath = `/${project_name}`
  return <Grid item container className={`${classes.header} ${className}`} >
    {
      link && projectPath != window.location.pathname ? <Link underline='none' className={ classes.link } href={projectPath}>{title}</Link> : title
    }
    {
      children
    }
  </Grid>
}

export default ProjectHeader;
