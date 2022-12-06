import React, {useContext} from 'react';
import Grid from '@material-ui/core/Grid';
import {makeStyles} from '@material-ui/core/styles';

import QueryFromPane from './query_from_pane';
import QuerySelectPane from './query_select_pane';
import QueryWherePane from './query_where_pane';

import {QueryGraphContext} from '../../contexts/query/query_graph_context';
import {QueryColumnContext} from '../../contexts/query/query_column_context';
import {QueryWhereContext} from '../../contexts/query/query_where_context';
import useUriQueryParams from '../../contexts/query/use_uri_query_params';

const useStyles = makeStyles((theme) => ({
  container: {
    width: '100%'
  },
  item: {
    width: '100%'
  }
}));

const QueryControls = ({
  syncQueryParams,
  queryControlButtons
}: {
  syncQueryParams: boolean;
  queryControlButtons: React.FC;
}) => {
  const classes = useStyles();
  const {
    state: {rootModel},
    setRootModel
  } = useContext(QueryGraphContext);
  const {state: columnState, setQueryColumns} = useContext(QueryColumnContext);
  const {state: whereState, setWhereState} = useContext(QueryWhereContext);

  useUriQueryParams({
    setRootModel,
    columnState,
    setQueryColumns,
    rootModel: rootModel || '',
    whereState,
    setWhereState,
    syncQueryParams
  });

  return (
    <Grid
      item
      container
      className={classes.container}
      justify='flex-start'
      alignItems='center'
      direction='column'
      xs={12}
    >
      <QueryFromPane queryControlButtons={queryControlButtons} />
      <QueryWherePane />
      <QuerySelectPane />
    </Grid>
  );
};

export default QueryControls;
