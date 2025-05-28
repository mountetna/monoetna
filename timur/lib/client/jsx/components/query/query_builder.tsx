import React, {useContext,useState,useRef} from 'react';

import {useReduxState} from 'etna-js/hooks/useReduxState';
import QueryControls from './query_controls';
import Grid from '@material-ui/core/Grid';
import Slide from '@material-ui/core/Slide';
import Tooltip from '@material-ui/core/Tooltip';
import QueryDisplayButtons from './query_display_buttons';
import QueryControlButtons from './query_control_buttons';
import QueryResults from './query_results';
import QueryString from './query_string';
import QueryOptions from './query_options';
import {QueryGraphContext} from '../../contexts/query/query_graph_context';

import useUriQueryParams from '../../contexts/query/use_uri_query_params';
import useQueryGraph from '../../contexts/query/use_query_graph';
import { useTheme, makeStyles } from '@material-ui/core/styles';

const useStyles = makeStyles((theme) => ({
  toolbar: {
    padding: '15px',
    borderBottom: '1px solid #ccc'
  },
  buttons: {
    borderRadius: '2px'
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
    overflowY: 'auto',
    flex: '1 1 calc(100% - max(700px, 40%))',
    flexWrap: 'nowrap',
    height: '100%'
  }
}));


const QueryBuilder = ({}) => {
  const {
    state: {graph},
    setGraph
  } = useContext(QueryGraphContext);
  const reduxState = useReduxState();

  useQueryGraph(reduxState, graph, setGraph);
  useUriQueryParams({});

  const classes = useStyles();

  const [ options, setOptions ] = useState(['controls']);

  const showQuery = options.includes('controls');
  const showRawQuery = options.includes('raw');

  const ref = useRef(null);

  const theme = useTheme();

  const transitionStyle = {
    transform: showQuery ? 'none' : 'translate(calc( -1 * max(700px,40%)))',
    flex: showQuery ? '1 1 max(700px, 40%)' : '0',
    transitionProperty: 'transform, flex',
    transitionDuration: `${theme.transitions.duration.standard}ms`,
    transitionTimingFunction: `${theme.transitions.easing.easeIn}`
  };

  if (!graph || !graph.initialized) return null;

  return (
      <Grid style={{ overflow: 'hidden', height: '100%' }} direction='column' container>
        <Grid item container direction='row' className={classes.toolbar} alignItems='center'>
          <QueryDisplayButtons
            options={options}
            setOptions={setOptions}/>
          <Grid item container direction='row' alignItems='center'
            justifyContent='flex-end'
            style={{ width: 'auto', flex: '1 1 auto' }}>
            <QueryControlButtons />
          </Grid>
        </Grid>
        <Grid item className={classes.query} container direction='row' ref={ref}>
            <Grid
              style={transitionStyle}
              item
              container
              className={classes.slide}
              justifyContent='flex-start'
              alignItems='center'
              direction='column'
            >
              { showRawQuery && <QueryString/> }
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
