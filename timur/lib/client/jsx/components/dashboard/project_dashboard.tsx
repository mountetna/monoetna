import React, {useState, useCallback, useEffect, useMemo} from 'react';

import {makeStyles} from '@material-ui/core/styles';
import Grid from '@material-ui/core/Grid';
import LinearProgress from '@material-ui/core/LinearProgress';
import Typography from '@material-ui/core/Typography';

import Button from '@material-ui/core/Button';
import List from '@material-ui/core/List';
import ListItem from '@material-ui/core/ListItem';
import ListItemText from '@material-ui/core/ListItemText';
import ListItemIcon from '@material-ui/core/ListItemIcon';
import HelpIcon from '@material-ui/icons/Help';
import StarIcon from '@material-ui/icons/Star';
import StarBorderIcon from '@material-ui/icons/StarBorder';
import StarHalfIcon from '@material-ui/icons/StarHalf';

import Card from '@material-ui/core/Card';
import CardContent from '@material-ui/core/CardContent';
import CardMedia from '@material-ui/core/CardMedia';

import {json_get, json_delete} from 'etna-js/utils/fetch';

const STAR = [
  StarBorderIcon,
  StarHalfIcon,
  StarIcon
];

const size = '35px';
const useStyles = makeStyles((theme) => ({
  app: {
    padding: '0px 15px',
    borderTop: '1px solid #eee',
    height: '90px',
  },
  appImage: {
    border: '1px solid green',
    borderRadius: size,
    backgroundSize: `${size} ${size}`,
    width: size,
    height: size,
  },
  appTitle: {
    fontSize: '0.9em',
    textTransform: 'capitalize',
    textAlign: 'center',
    color: theme.palette.secondary.main
  },
  list: {
    width: '100%',
    maxWidth: '36ch',
    backgroundColor: theme.palette.background.paper,
    padding: '0px'
  },
  list_item: {
    padding: '0px'
  },
  dashboard: {
    width: '600px',
    height: '600px'
  }
}));

const APPLICATIONS = {
}

type Info = {
  level: number;
  text: string;
}

const AppDashboard = ({app,title,action,children}:{
  app: string;
  title: string;
  action: string;
  children: any;
}) => {
  const classes = useStyles();

  return <Grid container alignContent='center' className={classes.app}>
    <Grid item direction='column' container alignContent='flex-start' alignItems='flex-start' justifyContent='space-evenly' style={{ width: '60px', padding: '0px', height: '100%', paddingBottom: '5px' }}>
      <Grid item
        style={{ width: '60px' }}><Typography className={classes.appTitle}>{title}</Typography></Grid>
      <Grid item style={{ paddingLeft: '10px' }}>
      <div
        className={classes.appImage}
        title={app}
        style={{
          backgroundImage: `url(/images/${app}.svg)`,
        }}
      />
      </Grid>
    </Grid>
    <Grid item container justifyContent='space-around' direction='column' style={{ width: 'auto', paddingLeft: '15px', flex: '1 1 auto' }}>
      <List className={classes.list}>
        {
          children
        }
      </List>
    </Grid>
    <Grid item container style={{width: '60px'}} alignContent='center' justifyContent='space-between' alignItems='center'>
      <HelpIcon style={{ color:'dodgerblue' }}/>
    </Grid>
  </Grid>
}

const STAR_STYLES = [
  { color: 'red' },
  { color: 'goldenrod' },
  { color: 'green' }
]

const AppInfo = ({sensor}:{sensor: Function})=> {
  const classes = useStyles();
  const [ info, setInfo ] = useState({ level: 0, text: undefined });

  useEffect( () => {
    sensor(setInfo);
  }, [] );

  const Star = STAR[info.level];

  return <ListItem className={classes.list_item}>
    <ListItemIcon>
      <Star style={STAR_STYLES[info.level]}/>
    </ListItemIcon>
    { info.text ?  <ListItemText secondary={info.text} /> : <div style={{width: '50px'}}><LinearProgress/></div> }
  </ListItem>
}

const JanusAppDashboard = ({project_name, mode}:{project_name: string, mode: any}) => {
  return <AppDashboard app='janus' title='access' action='Add users'>
    <AppInfo
      sensor={ (setInfo: Function) => json_get(`${CONFIG.janus_host}/api/admin/${project_name}/info`).then(
        ({project}) => {
          const summary = project.permissions.map( ({role}:{role: string}) => role ).reduce(
            (s: any,r: string) => {
              s[r] = (s[r] || 0) + 1;
              s.count = (s.count || 0) + 1;
              return s
            }, {}
          );
          const count = (c: number,r: string) => c == 1 ? `1 ${r}` : `${c} ${r}s`;
          const level = (summary.count < 3) ? (summary.count < 2 ? 0 : 1) : 2;
          const text = [ 'administrator', 'editor', 'viewer', 'guest' ].map(
            role => role in summary ? count(summary[role], role) : ''
          ).filter(_=>_).join(', ')
          setInfo({level, text})
        }
      ) }
    />
  </AppDashboard>
}

const TimurAppDashboard = ({project_name, mode}:{project_name: string, mode: any}) => {
  return <AppDashboard app='timur' title='modeling' action='Model'>
    <AppInfo sensor={ (setInfo: Function) => setInfo(
      { level: 1, text: '1 out of 3 models with attributes' }
    ) }/>
  </AppDashboard>
}

const GnomonAppDashboard = ({project_name, mode}:{project_name: string, mode: any}) => {
  return <AppDashboard app='gnomon' title='naming' action='Edit rules'>
    <AppInfo sensor={ (setInfo: Function) => setInfo(
      { level: 0, text: '1 of 12 models have identifier rules' }
    ) }/>
  </AppDashboard>
}

const MetisAppDashboard = ({project_name, mode}:{project_name: string, mode: any}) => {
  return <AppDashboard app='metis' title='files' action='Add files'>
    <AppInfo sensor={ (setInfo: Function) => setInfo(
      { level: 2, text: '2 buckets created' }
    ) }/>
    <AppInfo sensor={ (setInfo: Function) => setInfo(
      { level: 2, text: '572 files, 2.03 TB stored'}
    ) }/>
  </AppDashboard>
}

const PolyphemusAppDashboard = ({project_name, mode}:{project_name: string, mode: any}) => {
  return <AppDashboard app='polyphemus' title='linking' action='Link records'>
    <AppInfo sensor={ (setInfo: Function) => setInfo(
      { level: 2, text: '200 records created' }
    ) }/>
    <AppInfo sensor={ (setInfo: Function) => setInfo(
      { level: 2, text: '1 data loader, last run 2025-02-02' }
    ) }/>
  </AppDashboard>
}

const VulcanAppDashboard = ({project_name, mode}:{project_name: string, mode: any}) => {
  return <AppDashboard app='vulcan' title='analysis' action='Add workflows'>
    <AppInfo sensor={ (setInfo: Function) => setInfo(
      { level: 2, text: '1 workflow' }
    ) }/>
    <AppInfo sensor={ (setInfo: Function) => setInfo(
      { level: 2, text: '1 workspace, last run 2025-03-03' }
    ) }/>
  </AppDashboard>
}

const Dashboard = ({project_name}:{project_name: string}) => {
  const mode = 'administrator';
  const classes = useStyles();
  return <Card className={ classes.dashboard }>
    <CardContent>
    <Typography>Dashboard</Typography>
    </CardContent>
    <CardMedia>
      <JanusAppDashboard project_name={project_name} mode={mode}/>
      <TimurAppDashboard project_name={project_name} mode={mode}/>
      <GnomonAppDashboard project_name={project_name} mode={mode}/>
      <MetisAppDashboard project_name={project_name} mode={mode}/>
      <PolyphemusAppDashboard project_name={project_name} mode={mode}/>
      <VulcanAppDashboard project_name={project_name} mode={mode}/>
    </CardMedia>
  </Card>
}

export default Dashboard;
