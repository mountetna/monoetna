import React, {useContext, useMemo} from 'react';
import Grid from '@material-ui/core/Grid';
import {makeStyles} from '@material-ui/core/styles';

import {QueryGraphContext} from '../../contexts/query/query_graph_context';
import {QueryColumnContext} from '../../contexts/query/query_column_context';
import {QueryResultsContext} from '../../contexts/query/query_results_context';
import {userColumns} from '../../selectors/query_selector';
import QueryTable from './query_table';
import useTableEffects from './query_use_table_effects';

const useStyles = makeStyles((theme) => ({
  checkbox: {
    marginRight: '30px'
  },
  results: {
    height: '100%',
    overflowY: 'hidden'
  }
}));

const QueryResults = () => {
  const {
    state: {graph, rootModel}
  } = useContext(QueryGraphContext);
  const {
    state: {columns}
  } = useContext(QueryColumnContext);
  const {
    state: {
      pageSize,
      page,
      numRecords,
      queryString,
      data,
      expandMatrices,
      maxColumns
    },
    setPageSize,
    setPage
  } = useContext(QueryResultsContext);

  const classes = useStyles();

  const {columns: formattedColumns, rows} = useTableEffects({
    columns,
    data,
    graph,
    expandMatrices,
    maxColumns
  });

  function handlePageSizeChange(
    e: React.ChangeEvent<HTMLTextAreaElement | HTMLInputElement>
  ) {
    setPageSize(parseInt(e.target.value));
  }

  function handlePageChange(
    e: React.MouseEvent<HTMLButtonElement, MouseEvent> | null,
    newPage: number
  ) {
    setPage(newPage);
  }

  const userColumnsStr = useMemo(() => {
    return JSON.stringify(userColumns(columns));
  }, [columns]);

  if (!rootModel) return null;

  return (
    <Grid xs={12} item container direction='column' className={classes.results}>
      <QueryTable
        maxColumns={maxColumns}
        columns={formattedColumns}
        rows={rows}
        pageSize={pageSize}
        numRecords={numRecords}
        page={page}
        expandMatrices={expandMatrices}
        handlePageChange={handlePageChange}
        handlePageSizeChange={handlePageSizeChange}
      />
    </Grid>
  );
};

export default QueryResults;
