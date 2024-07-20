import React, {useContext,useState,useRef} from 'react';

import {useReduxState} from 'etna-js/hooks/useReduxState';
import QueryControls from './query_controls';
import Grid from '@mui/material/Grid';
import Slide from '@mui/material/Slide';
import ToggleButton from '@mui/material/ToggleButton';
import MenuIcon from '@mui/icons-material/Menu';
import QueryControlButtons from './query_control_buttons';
import QueryResults from './query_results';
import {QueryGraphContext} from '../../contexts/query/query_graph_context';

import useQueryGraph from '../../contexts/query/use_query_graph';
import { makeStyles } from '@mui/styles';

const useStyles = makeStyles((theme) => ({
  buttons: {
    padding: '15px',
    borderBottom: '1px solid #ccc'
  },
  slide: {
    overflowY: 'scroll',
    height: '100%',
    flexWrap: 'nowrap',
    borderRight: '1px solid #ccc'
  },
  container: {
    fontFamily: 'Baskervville',
    fontWeight: 100,
    fontSize: '1.2em',
    width: 'auto',
    height: '100%',
    overflowY: 'scroll'
  },
  query: {
    flex: '0 1 calc(100% - 79px)',
    overflow: 'hidden'
  },
  results: {
    boxShadow: '1px 1px 5px -2px #ccc inset',
    overflowX: 'auto',
    flex: '1 1 calc(100% - max(700px, 40%))',
    flexWrap: 'nowrap'
  }
}));


const QueryBuilder = ({}) => {
  const {
    state: {graph},
    setGraph
  } = useContext(QueryGraphContext);
  const reduxState = useReduxState();

  useQueryGraph(reduxState, graph, setGraph);

  if (!graph || !graph.initialized) return null;

  const classes = useStyles();

  const [ showQuery, setShowQuery ] = useState(true);

  const ref = useRef(null);

  const transitionStyle = theme => ({
    transform: showQuery ? 'none' : 'translate(calc( -1 * max(700px,40%)))',
    flex: showQuery ? '1 1 max(700px, 40%)' : '0',
    transitionProperty: "transform, flex",
    transitionDuration: `${theme.transitions.duration.standard}ms`,
    trasitionTimingFunction: `${theme.transitions.easing.easeIn}`
  });

  return (
      <Grid style={{ overflow: 'hidden', height: '100%' }} direction='column' container>
        <Grid item container direction='row' className={classes.buttons} alignItems='center'>
          <Grid item><ToggleButton selected={showQuery} onChange={ () => setShowQuery(!showQuery) } color='secondary'><MenuIcon/></ToggleButton></Grid>
          <Grid item container direction='row' alignItems='center'
            justifyContent='flex-end'
            style={{ width: 'auto', flex: '1 1 auto' }}>
            <QueryControlButtons />
          </Grid>
        </Grid>
        <Grid item className={classes.query} container direction='row' ref={ref}>
            <Grid
              sx={transitionStyle}
              item
              container
              className={classes.slide}
              justifyContent='flex-start'
              alignItems='center'
              direction='column'
            >
              <QueryControls />
            </Grid>
          <Grid direction='column' container className={classes.results}>
            <QueryResults />
          </Grid>
        </Grid>
      </Grid>
  );
};

export default QueryBuilder;
