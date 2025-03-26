import React, {useState, useCallback, useEffect, useMemo} from 'react';

import ProjectDashboard from './project_dashboard';
import ModelMapGraphic from '../model_map/model_map_graphic';

import {makeStyles} from '@material-ui/core/styles';
import Grid from '@material-ui/core/Grid';

const useStyles = makeStyles((theme) => ({
}));

const Dashboard = ({project_name}) => {
  return <Grid container style={{ height: 'calc(100vh - 150px)', width: '100vw', flexWrap: 'nowrap' }}>
    <Grid item container style={{ overflow: 'auto', width: 'calc(100vw - 600px)' }}>
      <ProjectDashboard project_name={project_name}/>
    </Grid>
    <Grid item container style={{ overflow: 'auto', width: '600px' }}>
      <ModelMapGraphic
        width={600}
        height={600}
        selected_models={[]}
      />
    </Grid>
  </Grid>
}

export default Dashboard;
