import React, {useState, createContext, useContext, useCallback} from 'react';

import {QueryResponse, SavedQuery} from './query_types';
import {QueryGraphContext} from '../../contexts/query/query_graph_context';
import {QueryColumnContext} from '../../contexts/query/query_column_context';
import {QueryWhereContext} from '../../contexts/query/query_where_context';
import {decodeCompressedParams, unpackParams} from '../../utils/query_uri_params';
import {QueryBuilder} from '../../utils/query_builder';

export const defaultQueryResultsParams = {
  expandMatrices: true,
  flattenQuery: true,
  showDisconnected: false,
  page: 0,
  pageSize: 10,
  data: {} as QueryResponse,
  numRecords: 0,
  queryString: '',
  savedQueries: [] as SavedQuery[],
  maxColumns: 10
};

const defaultQueryResultsState = {
  ...defaultQueryResultsParams
};

export type QueryResultsState = Readonly<typeof defaultQueryResultsState>;

export const defaultQueryResultsContext = {
  state: defaultQueryResultsState as QueryResultsState,
  setExpandMatrices: (expandMatrices: boolean) => {},
  setFlattenQuery: (flattenQuery: boolean) => {},
  setShowDisconnected: (showDisconnected: boolean) => {},
  setPage: (page: number) => {},
  setPageSize: (pageSize: number) => {},
  setDataAndNumRecords: (data: QueryResponse, numRecords: number) => {},
  setQueryString: (queryString: string) => {},
  setSavedQueries: (queries: SavedQuery[]) => {},
  setResultsState: (newState: QueryResultsState) => {},
  setQueryStateFromString: async (search: string) => {}
};

export type QueryResultsContextData = typeof defaultQueryResultsContext;

export const QueryResultsContext = createContext(defaultQueryResultsContext);
export type QueryResultsContext = typeof QueryResultsContext;
export type ProviderProps = {params?: {}; children: any};

export const QueryResultsProvider = (
  props: ProviderProps & Partial<QueryResultsContextData>
) => {
  const [state, setState] = useState(
    props.state || defaultQueryResultsContext.state
  );

  const {
    state: {rootModel='', graph},
    setRootModel
  } = useContext(QueryGraphContext);
  const {state: columnState, setQueryColumns} = useContext(QueryColumnContext);
  const {state: whereState, setWhereState} = useContext(QueryWhereContext);

  const setExpandMatrices = useCallback(
    (expandMatrices: boolean) => {
      setState({
        ...state,
        expandMatrices
      });
    },
    [state]
  );

  const setFlattenQuery = useCallback(
    (flattenQuery: boolean) => {
      setState({
        ...state,
        flattenQuery
      });
    },
    [state]
  );

  const setShowDisconnected = useCallback(
    (showDisconnected: boolean) => {
      setState({
        ...state,
        showDisconnected
      });
    },
    [state]
  );

  const setPage = useCallback(
    (page: number) => {
      setState({
        ...state,
        page
      });
    },
    [state]
  );

  const setPageSize = useCallback(
    (pageSize: number) => {
      setState({
        ...state,
        pageSize
      });
    },
    [state]
  );

  const setDataAndNumRecords = useCallback(
    (data: QueryResponse, numRecords: number) => {
      setState({
        ...state,
        data,
        numRecords
      });
    },
    [state]
  );

  const setQueryString = useCallback(
    (queryString: string) => {
      setState({
        ...state,
        queryString
      });
    },
    [state]
  );

  const unpackQuery = async (query: SavedQuery) => {
    if (query.unpackedQuery) return query;

    const queryState = await decodeCompressedParams(query.query);

    const builder = new QueryBuilder(graph);
    builder.addRootModel(queryState.rootModel);
    builder.addColumns(queryState.columns);
    builder.addRecordFilters(queryState.recordFilters);
    builder.setFlatten(queryState.flattenQuery);
    builder.setOrRecordFilterIndices(queryState.orRecordFilterIndices);

    const unpackedQuery = JSON.stringify(builder.query(), null, 2);

    return {
      ...query,
      unpackedQuery 
    }
  };

  const setSavedQueries = useCallback(
    async (savedQueries: SavedQuery[]) => {
      setState({
        ...state,
        savedQueries: await Promise.all(savedQueries.map(unpackQuery))
      });
    },
    [state]
  );

  const setResultsState = useCallback((newState: QueryResultsState) => {
    setState({
      ...newState
    });
  }, []);

  const setQueryStateFromString = async (search: string) => {
    let params = await unpackParams(search);

    setQueryColumns(params.columns);
    setRootModel(params.rootModel);
    setWhereState({
      recordFilters: params.recordFilters,
      orRecordFilterIndices: params.orRecordFilterIndices
    });
  }

  return (
    <QueryResultsContext.Provider
      value={{
        state,
        setQueryStateFromString,
        setExpandMatrices,
        setFlattenQuery,
        setShowDisconnected,
        setPage,
        setPageSize,
        setDataAndNumRecords,
        setQueryString,
        setSavedQueries,
        setResultsState
      }}
    >
      {props.children}
    </QueryResultsContext.Provider>
  );
};
