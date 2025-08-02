import React, {useState, useCallback, useEffect, useMemo} from 'react';

import ProjectDashboard from './project_dashboard';
import DashboardMap from './dashboard_map';
import ModelMapGraphic from '../model_map/model_map_graphic';
import {DashboardProvider} from '../../contexts/dashboard_context';

import { useTheme, makeStyles } from '@material-ui/core/styles';
import IconButton from '@material-ui/core/IconButton';
import BuildIcon from '@material-ui/icons/Build';
import MapIcon from '@material-ui/icons/Map';
import ChevronRightIcon from '@material-ui/icons/ChevronRight';
import CloseIcon from '@material-ui/icons/Close';
import Grid from '@material-ui/core/Grid';
import Drawer from '@material-ui/core/Drawer';
import Tabs from '@material-ui/core/Tabs';
import Tab from '@material-ui/core/Tab';
import Tooltip from '@material-ui/core/Tooltip';

const useStyles = makeStyles((theme) => ({
  tab: {
    minWidth: '60px'
  },
  controls: {
    width: '60px',
    boxShadow: '0 0 5px #888',
    clipPath: 'inset(0 0 -10px -10px)',
    zIndex: 1,
    background: '#fbf9f6'
  },
  tabs: {
    background: 'goldenrod'
  },
  tab_panel: {
    background: 'white',
    boxShadow: '0 0 5px #888',
    clipPath: 'inset(0 -10px -10px -10px)',
  },
  topshadow: {
    boxShadow: '10px 0 5px #888',
    clipPath: 'inset(-10px 0 0 0)',
    paddingLeft: '10px'
  },
  slide: {
    overflow: 'visible',
    height: '650px',
    width: '672px',
    position: 'fixed',
    top: '70px',
    right: '-1px',
    clipPath: 'inset(-10px 0 0 -10px)'
  },
}));

type ModelState = [ modelName: string|null, modelEl: HTMLElement|null ];

const Dashboard = ({project_name, model_name}:{project_name: string, model_name?: string}) => {
  const [ showDrawer, setShowDrawer ] = useState(model_name == 'project');
  const [ tab, setTab ] = useState(0);
  const [ modelState, setModelState ] = useState<ModelState>([null, null]);

  const classes = useStyles();
  const theme = useTheme();

  const transitionStyle = {
    transform: showDrawer ? 'none' : 'translate(604px)',
    transitionProperty: 'transform',
    transitionDuration: `${theme.transitions.duration.standard}ms`,
    transitionTimingFunction: `${theme.transitions.easing.easeIn}`
  };

  return <DashboardProvider>
    <Grid
        style={transitionStyle}
        item
        container
        className={classes.slide}
        justifyContent='flex-start'
      >
        <Grid container className={classes.topshadow} alignItems='flex-start' >
          <Grid container direction='column' className={ classes.controls }>
            <Tabs value={tab} onChange={ (e, tab) => {
              setTab(tab); setShowDrawer(true);
            } } orientation='vertical'
            className={ classes.tabs }>
              <Tab className={ classes.tab } disableRipple icon={
                <Tooltip title="dashboard"><BuildIcon/></Tooltip>
              } />
              <Tab className={ classes.tab } disableRipple icon={
                <Tooltip title="map"><MapIcon/></Tooltip>
              } />
            </Tabs>
            {
              showDrawer && <IconButton onClick={ () => { setShowDrawer(false); setModelState([null,null]); } }><CloseIcon/></IconButton>
            }
          </Grid>
          <Grid className={classes.tab_panel}>
          {
            tab == 0 && <ProjectDashboard project_name={project_name}/>
          }
          {
            tab == 1 && <DashboardMap modelState={modelState} setModelState={setModelState}/>
          }
        </Grid>
      </Grid>
    </Grid>
  </DashboardProvider>
}

export default Dashboard;
