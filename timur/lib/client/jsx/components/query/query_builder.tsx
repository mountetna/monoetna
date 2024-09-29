import React, {useContext,useState,useRef} from 'react';

import {useReduxState} from 'etna-js/hooks/useReduxState';
import QueryControls from './query_controls';
import Grid from '@material-ui/core/Grid';
import Slide from '@material-ui/core/Slide';
import Tooltip from '@material-ui/core/Tooltip';
import IconButton from '@material-ui/core/IconButton';
import ButtonGroup from '@material-ui/core/ButtonGroup';
import ListIcon from '@material-ui/icons/List';
import CodeIcon from '@material-ui/icons/Code';
import TuneIcon from '@material-ui/icons/Tune';
import QueryControlButtons from './query_control_buttons';
import QueryResults from './query_results';
import QueryString from './query_string';
import QueryOptions from './query_options';
import {QueryGraphContext} from '../../contexts/query/query_graph_context';

import useQueryGraph from '../../contexts/query/use_query_graph';
import { useTheme, makeStyles } from '@material-ui/core/styles';

const useStyles = makeStyles((theme) => ({
  toolbar: {
    padding: '15px',
    borderBottom: '1px solid #ccc'
  },
  button: {
    '&& :not(:last-child)': {
      borderRight: '1px solid #ccc'
    }
  },
  buttons: {
    border: '1px solid #ccc',
    borderRadius: '2px',
    height: '30px'
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
  const [ showRawQuery, setShowRawQuery ] = useState(false);
  const [ showQueryOptions, setShowQueryOptions ] = useState(false);

  const ref = useRef(null);

  const theme = useTheme();

  const transitionStyle = {
    transform: showQuery ? 'none' : 'translate(calc( -1 * max(700px,40%)))',
    flex: showQuery ? '1 1 max(700px, 40%)' : '0',
    transitionProperty: "transform, flex",
    transitionDuration: `${theme.transitions.duration.standard}ms`,
    transitionTimingFunction: `${theme.transitions.easing.easeIn}`
  };

  return (
      <Grid style={{ overflow: 'hidden', height: '100%' }} direction='column' container>
        <Grid item container direction='row' className={classes.toolbar} alignItems='center'>
          <ButtonGroup className={classes.buttons} size="small">
            <Tooltip title="Query controls">
              <IconButton className={classes.button} disableRipple={true} color={ showQuery ? 'secondary' : 'primary' } onClick={ () => setShowQuery(!showQuery) }>
                <ListIcon/>
              </IconButton>
            </Tooltip>
            { showQuery && <Tooltip title="Raw query">
              <IconButton className={classes.button} disableRipple={true} color={showRawQuery ? 'secondary' : 'primary' } onClick={ () => setShowRawQuery(!showRawQuery) } >
                  <CodeIcon/>
                </IconButton>
              </Tooltip>
            }
            { showQuery && <Tooltip title="Query Options">
                <IconButton className={classes.button} disableRipple={true} color={showQueryOptions ? 'secondary' : 'primary' } onClick={ () => setShowQueryOptions(!showQueryOptions) }>
                  <TuneIcon/>
                </IconButton>
              </Tooltip>
            }
          </ButtonGroup>
          <QueryOptions open={showQueryOptions} onClose={ () => setShowQueryOptions(false) }/>
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
