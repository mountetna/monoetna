import React, {useCallback, useContext} from 'react';

import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
import IconButton from '@material-ui/core/IconButton';
import {makeStyles} from '@material-ui/core/styles';

import {QueryGraphContext} from '../../contexts/query/query_graph_context';
import {columnsForModel, QueryColumnContext} from '../../contexts/query/query_column_context';
import {QueryWhereContext} from '../../contexts/query/query_where_context';
import {QueryResultsContext} from '../../contexts/query/query_results_context';
import {EmptyQueryResponse} from '../../contexts/query/query_types';
import QueryModelSelector from './query_model_selector';
import QueryClause from './query_clause';
import QueryChevron from './query_chevron';
import {Attribute, Model} from 'etna-js/models/magma-model';

const useStyles = makeStyles((theme) => ({
  title: {
    fontWeight: 'bold'
  },
  summary: {
    fontStyle: 'italic',
    paddingLeft: '10px',
  }
}));
const QueryRowPane = () => {
  const {
    state: {graph, rootModel},
    setRootModel
  } = useContext(QueryGraphContext);
  const {setQueryColumns} = useContext(QueryColumnContext);
  const {resetWhereState} = useContext(QueryWhereContext);
  const {setDataAndNumRecords} = useContext(QueryResultsContext);

  const onRootModelSelect = useCallback(
    (modelName: string) => {
      let model: Model = graph.models.model(modelName) as Model;
      setRootModel(modelName);
      resetWhereState();
      setQueryColumns(columnsForModel(model, graph.models));
      setDataAndNumRecords(EmptyQueryResponse, 0);
    },
    [graph, setRootModel, setQueryColumns, setDataAndNumRecords]
  );

  const classes = useStyles();

  return (
    <QueryClause title=''>
      <Grid item container>
        <Grid item container xs={8} alignItems='center'>
          <QueryChevron disabled/>
          <Typography className={classes.title}>Rows:</Typography>
          <Typography className={classes.summary}>from&nbsp;</Typography>
          <QueryModelSelector
            setModel={onRootModelSelect}
            modelNames={graph.models.modelNames}
            modelName={rootModel || ''}
          />
        </Grid>
      </Grid>
    </QueryClause>
  );
};

export default QueryRowPane;
