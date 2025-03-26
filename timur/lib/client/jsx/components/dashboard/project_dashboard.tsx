import React, {useState, useCallback, useEffect, useMemo} from 'react';

import {makeStyles} from '@material-ui/core/styles';
import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';

import List from '@material-ui/core/List';
import ListItem from '@material-ui/core/ListItem';
import ListItemIcon from '@material-ui/core/ListItemIcon';
import ListItemText from '@material-ui/core/ListItemText';
import Divider from '@material-ui/core/Divider';

const size = '35px';
const useStyles = makeStyles((theme) => ({
  app: {
    padding: '0px 15px',
    borderBottom: '1px solid #eee'
  },
  appImage: {
    border: '2px solid #066306',
    borderRadius: size,
    backgroundSize: `${size} ${size}`,
    width: size,
    height: size,
  },
  appTitle: {
    fontSize: '0.9em',
    textTransform: 'capitalize',
    color: theme.palette.secondary.main
  },
  list: {
    width: '100%',
    maxWidth: '36ch',
    backgroundColor: theme.palette.background.paper
  },
}));

const APPLICATIONS = {
}

const AppDashboard = ({app,title,items}) => {
  const classes = useStyles();
  return <Grid container alignContent='center' className={classes.app}>
    <Grid item direction='column' container alignContent='center' alignItems='center' justifyContent='flex-start' style={{ width: '60px', padding: '10px' }}>
      <Grid item>
      <div
        className={classes.appImage}
        title={app}
        style={{
          backgroundImage: `url(/images/${app}.svg)`,
        }}
      />
      </Grid>
      <Grid item><Typography className={classes.appTitle}>{title}</Typography></Grid>
    </Grid>
    <Grid container alignContent='center' direction='column' style={{ width: 'auto', paddingLeft: '25px' }}>
      <List className={classes.list}>
        {
          items.map( (item,i) =>
            <ListItem key={i}>
              <ListItemText primary={item} />
            </ListItem>
          )
        }
      </List>
    </Grid>
  </Grid>
}

const JanusAppDashboard = ({project_name, mode}) => {
  return <AppDashboard app='janus' title='access' items={[ '2 administrators' ]}/>
}

const TimurAppDashboard = ({project_name, mode}) => {
  return <AppDashboard app='timur' title='modeling' items={[ '2 models with attributes', '1 empty model' ]}/>
}

const GnomonAppDashboard = ({project_name, mode}) => {
  return <AppDashboard app='gnomon' title='naming' items={[ '1 of 12 models have identifier rules' ]}/>
}

const MetisAppDashboard = ({project_name, mode}) => {
  return <AppDashboard app='metis' title='files' items={[ '572 files', '2033 GB' ]}/>
}

const PolyphemusAppDashboard = ({project_name, mode}) => {
  return <AppDashboard app='polyphemus' title='linking' items={[ '200 records created', '1 data loader', 'Last loader ran 2025-02-02' ]}/>
}

const VulcanAppDashboard = ({project_name, mode}) => {
  return <AppDashboard app='vulcan' title='analysis' items={[ '1 analysis workspace' ]}/>
}

const Dashboard = ({project_name}) => {
  const items = [];
  const mode = 'administrator';
  return <Grid container alignContent='flex-start'>
    <JanusAppDashboard project_name={project_name} mode={mode}/>
    <TimurAppDashboard project_name={project_name} mode={mode}/>
    <GnomonAppDashboard project_name={project_name} mode={mode}/>
    <MetisAppDashboard project_name={project_name} mode={mode}/>
    <PolyphemusAppDashboard project_name={project_name} mode={mode}/>
    <VulcanAppDashboard project_name={project_name} mode={mode}/>
  </Grid>
}

export default Dashboard;
