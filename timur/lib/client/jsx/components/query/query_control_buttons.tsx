import React, {
  useMemo,
  useContext,
  useEffect,
  useState,
  useCallback
} from 'react';
import * as _ from 'lodash';
import Button from '@material-ui/core/Button';
import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
import ReplayIcon from '@material-ui/icons/Replay';
import ShareIcon from '@material-ui/icons/Share';

import {makeStyles} from '@material-ui/core/styles';

import {useFeatureFlag} from 'etna-js/hooks/useFeatureFlag';
import {useModal} from 'etna-js/components/ModalDialogContainer';
import {copyText} from 'etna-js/utils/copy';
import {QueryGraphContext} from '../../contexts/query/query_graph_context';
import {
  QueryColumnContext,
  defaultQueryColumnParams
} from '../../contexts/query/query_column_context';
import {
  QueryWhereContext,
  defaultQueryWhereParams
} from '../../contexts/query/query_where_context';
import {
  QueryResultsContext,
  defaultQueryResultsParams
} from '../../contexts/query/query_results_context';
import {QueryBuilder} from '../../utils/query_builder';
import useTableEffects from './query_use_table_effects';
import useResultsActions from './query_use_results_actions';
import QueryTsvOptionsModal from './query_tsv_options_modal';
import SaveQueryModal from './save_query_modal';

const useStyles = makeStyles((theme) => ({
  button: {
    marginLeft: '5px'
  }
}));

const QueryControlButtons = () => {
  const {
    state: {graph, rootModel},
    setRootModel
  } = useContext(QueryGraphContext);
  const {
    state: {columns},
    setQueryColumns
  } = useContext(QueryColumnContext);
  const {
    state: {recordFilters, orRecordFilterIndices},
    setWhereState
  } = useContext(QueryWhereContext);
  const {
    state: {
      pageSize,
      page,
      data,
      expandMatrices,
      flattenQuery,
      showDisconnected,
      queryString,
      maxColumns
    },
    setQueryString,
    setPage,
    setDataAndNumRecords,
    setResultsState
  } = useContext(QueryResultsContext);

  const [lastColumns, setLastColumns] = useState(
    defaultQueryColumnParams.columns
  );
  const [lastFilters, setLastFilters] = useState(
    defaultQueryWhereParams.recordFilters
  );
  const [lastOrFilterIndices, setLastOrFilterIndices] = useState(
    defaultQueryWhereParams.orRecordFilterIndices
  );
  const [lastExpandMatrices, setLastExpandMatrices] = useState(expandMatrices);
  const [lastShowDisconnected, setLastShowDisconnected] = useState(showDisconnected);
  const [lastFlattenQuery, setLastFlattenQuery] = useState(flattenQuery);
  const [lastPage, setLastPage] = useState(page);
  const [lastPageSize, setLastPageSize] = useState(pageSize);

  const {openModal} = useModal();

  const classes = useStyles();

  const builder = useMemo(() => {
    if (rootModel && graph && graph.initialized) {
      let builder = new QueryBuilder(graph);

      builder.addRootModel(rootModel);
      builder.addColumns(columns);
      builder.addRecordFilters(recordFilters);
      builder.setFlatten(flattenQuery);
      builder.setOrRecordFilterIndices(orRecordFilterIndices);

      return builder;
    }

    return null;
  }, [
    columns,
    recordFilters,
    graph,
    orRecordFilterIndices,
    flattenQuery,
    rootModel
  ]);

  const query = useMemo(() => {
    if (!builder) return '';

    return builder.query();
  }, [builder]);

  const count = useMemo(() => {
    if (!builder) return '';

    return builder.count();
  }, [builder]);

  const {columns: formattedColumns} = useTableEffects({
    columns,
    data,
    graph,
    expandMatrices,
    maxColumns
  });

  const {runQuery, downloadData} = useResultsActions({
    countQuery: count,
    query,
    page,
    pageSize,
    columns,
    expandMatrices,
    showDisconnected,
    setDataAndNumRecords
  });

  useEffect(() => {
    // At some point, we can probably cache data and only
    //   fetch when needed?
    if (
      lastPage !== page ||
      lastPageSize !== pageSize
    ) {
      runQuery();
      setLastPage(page);
      setLastPageSize(pageSize);
    }
  }, [
    page,
    pageSize,
    lastPage,
    lastPageSize,
    runQuery
  ]);

  useEffect(() => {
    if (JSON.stringify(query) !== queryString)
      {setQueryString(JSON.stringify(query));}
  }, [query, setQueryString, queryString]);

  function resetQuery() {
    setRootModel(null);
    setResultsState(defaultQueryResultsParams);
    setQueryColumns(defaultQueryColumnParams.columns);
    setWhereState(defaultQueryWhereParams);
  }

  function copyLink() {
    copyText(window.location.href);
  }

  const [ showSaveQuery, setShowSaveQuery ] = useState(false);

  const disableQueryBtn = useMemo(() => {
    return !rootModel || (
      _.isEqual(columns, lastColumns) &&
      _.isEqual(recordFilters, lastFilters) &&
      _.isEqual(orRecordFilterIndices, lastOrFilterIndices) &&
      lastShowDisconnected === showDisconnected &&
      lastExpandMatrices === expandMatrices &&
      lastFlattenQuery === flattenQuery
    );
  }, [
    columns,
    lastColumns,
    recordFilters,
    lastFilters,
    orRecordFilterIndices,
    lastOrFilterIndices,
    showDisconnected,
    lastShowDisconnected,
    flattenQuery,
    lastFlattenQuery,
    expandMatrices,
    lastExpandMatrices
  ]);

  const handleRunQuery = useCallback(() => {
    runQuery(0);
    setPage(0);
    setLastPage(0);
    setLastColumns(columns);
    setLastFilters(recordFilters);
    setLastOrFilterIndices(orRecordFilterIndices);
    setLastExpandMatrices(expandMatrices);
    setLastShowDisconnected(showDisconnected);
    setLastFlattenQuery(flattenQuery);
  }, [runQuery, columns, recordFilters, orRecordFilterIndices,expandMatrices,showDisconnected,flattenQuery]);

  return (
    <>
      <Grid item>
        <Button
          className={classes.button}
          onClick={resetQuery}
          startIcon={<ReplayIcon />}
        >
          Reset Query
        </Button>
        <Button
          className={classes.button}
          color='secondary'
          onClick={handleRunQuery}
          disabled={disableQueryBtn}
        >
          Query
        </Button>
        <Button
          className={classes.button}
          onClick={() => {
            openModal(<QueryTsvOptionsModal onDownload={downloadData} />);
          }}
        >
          {'\u21af TSV'}
        </Button>
        <Button
          className={classes.button}
          onClick={copyLink}
          startIcon={<ShareIcon />}
        >
          Copy Link
        </Button>
        <Button
          className={classes.button}
          color='secondary'
          onClick={ () => setShowSaveQuery(true) }
        >
          Save Query
        </Button>
        <SaveQueryModal
          open={ showSaveQuery }
          onClose={ () => setShowSaveQuery(false) }/>
      </Grid>
    </>
  );
};

export default QueryControlButtons;
