import React, {useState, useCallback, useEffect, useMemo} from 'react';

import {makeStyles} from '@material-ui/core/styles';
import Grid from '@material-ui/core/Grid';
import LinearProgress from '@material-ui/core/LinearProgress';
import Typography from '@material-ui/core/Typography';

import Button from '@material-ui/core/Button';
import IconButton from '@material-ui/core/IconButton';
import Tooltip from '@material-ui/core/Tooltip';
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

import {useReduxState} from 'etna-js/hooks/useReduxState';
import {selectUser} from 'etna-js/selectors/user-selector';

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
  link: {
    color: 'darkgoldenrod',
    cursor: 'pointer'
  },
  dashboard: {
    width: '600px',
    height: '100%'
  }
}));

const APPLICATIONS = {
}

type Info = {
  level: number;
  text: string;
}

const AppDashboard = ({app,title,action,help,helpLink,children}:{
  app: string;
  title: string;
  action: string;
  children: any;
}) => {
  const classes = useStyles();

  const user = useReduxState((state) => selectUser(state));
  const userRole = 'viewer';//user.permissions[CONFIG.project_name].role;

  const shownChildren = React.Children.map(
    children, child => ROLES[child.props.role] > ROLES[userRole] ? null : child
  ).filter(_=>_)

  if (!shownChildren.length) return null;

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
          shownChildren
        }
      </List>
    </Grid>
    <Grid item container style={{width: '60px'}} alignContent='center' justifyContent='space-between' alignItems='center'>
      <Tooltip title={help}><IconButton disableRipple href={helpLink}><HelpIcon style={{ color:'dodgerblue' }}/></IconButton></Tooltip>
    </Grid>
  </Grid>
}

const STAR_STYLES = [
  { color: 'red' },
  { color: 'goldenrod' },
  { color: 'green' }
]

const ROLES = {
  guest: 0,
  viewer: 1,
  editor: 2,
  administrator: 3
};

const AppInfo = ({sensor,action,role,action_role=role,action_link}:{sensor: Function, action: string, role: string})=> {
  const classes = useStyles();
  const [ info, setInfo ] = useState({ level: 0, text: undefined });

  const user = useReduxState((state) => selectUser(state));
  const userRole =  'viewer';//user.permissions[CONFIG.project_name].role;

  useEffect( () => {
    sensor(setInfo);
  }, [] );

  const Star = STAR[info.level];

  return <ListItem className={classes.list_item}>
    <ListItemIcon>
      <Star style={STAR_STYLES[info.level]}/>
    </ListItemIcon>
    { info.text
      ? <ListItemText secondary={
          ROLES[action_role] > ROLES[userRole]
            ? info.text
            : <a className={classes.link} title={action} href={action_link}>{info.text}</a>
        } />
      : <div style={{width: '50px'}}><LinearProgress/></div> }
  </ListItem>
}

const JanusAppDashboard = ({project_name}:{project_name: string}) => {
  return <AppDashboard app='janus' title='access' help='Guide to Managing User Access' helpLink='https://mountetna.github.io/access.html'>
    <AppInfo action='Add users'
      role='administrator'
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

const TimurAppDashboard = ({project_name}:{project_name: string}) => {
  return <AppDashboard app='timur' title='modeling' help='Guide to Modeling' helpLink='https://mountetna.github.io/modeling.html'>
    <AppInfo
      role='editor'
      action_role='administrator'
      action_link={ `/${project_name}/map` }
      action='Edit models' sensor={ (setInfo: Function) => setInfo(
      { level: 1, text: '1 out of 3 models with attributes' }
    ) }/>
  </AppDashboard>
}

const GnomonAppDashboard = ({project_name}:{project_name: string}) => {
  return <AppDashboard app='gnomon' title='naming' help='Guide to Naming' helpLink='https://mountetna.github.io/naming.html'>
    <AppInfo
      role='editor'
      action_link={ `https://${CONFIG.gnomon_host}/${project_name}/rules` }
      action='Edit rules' sensor={ (setInfo: Function) => setInfo(
      { level: 0, text: '1 of 12 models have identifier rules' }
    ) }/>
  </AppDashboard>
}

const MetisAppDashboard = ({project_name}:{project_name: string}) => {
  return <AppDashboard app='metis' title='files' help='Guide to Ingesting Files' helpLink='https://mountetna.github.io/ingestion.html'>
    <AppInfo
      role='editor'
      action_role='administrator'
      action_link={`https://${CONFIG.metis_host}/${project_name}/`}
      action='Create buckets' sensor={ (setInfo: Function) => setInfo(
      { level: 2, text: '2 buckets created' }
    ) }/>
    <AppInfo
      role='viewer'
      action_role='editor'
      action_link={`https://${CONFIG.metis_host}/${project_name}/`}
      action='Add files' sensor={ (setInfo: Function) => setInfo(
      { level: 2, text: '572 files, 2.03 TB stored'}
    ) }/>
  </AppDashboard>
}

const PolyphemusAppDashboard = ({project_name}:{project_name: string}) => {
  return <AppDashboard app='polyphemus' title='linking' help='Guide to Linking Records' helpLink='https://mountetna.github.io/linking.html'>
    <AppInfo
      role='viewer'
      action_role='editor'
      action_link={`https://${CONFIG.polyphemus_host}/${project_name}`}
      action='Link records' sensor={ (setInfo: Function) => setInfo(
      { level: 2, text: '200 records created' }
    ) }/>
    <AppInfo
      role='editor'
      action_link={`https://${CONFIG.polyphemus_host}/${project_name}`}
      action='Create loaders' sensor={ (setInfo: Function) => setInfo(
      { level: 2, text: '1 data loader, last run 2025-02-02' }
    ) }/>
  </AppDashboard>
}

const VulcanAppDashboard = ({project_name}:{project_name: string}) => {
  return <AppDashboard app='vulcan' title='analysis' help='Guide to Analysis Workflows' helpLink='https://mountetna.github.io/analysis.html'>
    <AppInfo
      action_role='administrator'
      action_link={`https://${CONFIG.vulcan_host}/${project_name}`}
      role='viewer' action='Add workflows' sensor={ (setInfo: Function) => setInfo(
      { level: 2, text: '1 workflow' }
    ) }/>
    <AppInfo
      role='viewer'
      action_link={`https://${CONFIG.vulcan_host}/${project_name}`}
      action='Run workflows' sensor={ (setInfo: Function) => setInfo(
      { level: 2, text: '1 workspace, last run 2025-03-03' }
    ) }/>
  </AppDashboard>
}

const Dashboard = ({project_name}:{project_name: string}) => {
  const classes = useStyles();
  return <Card className={ classes.dashboard }>
    <CardContent>
    <Typography>Dashboard</Typography>
    </CardContent>
    <CardMedia>
      <JanusAppDashboard project_name={project_name}/>
      <TimurAppDashboard project_name={project_name}/>
      <GnomonAppDashboard project_name={project_name}/>
      <MetisAppDashboard project_name={project_name}/>
      <PolyphemusAppDashboard project_name={project_name}/>
      <VulcanAppDashboard project_name={project_name}/>
    </CardMedia>
  </Card>
}

export default Dashboard;
