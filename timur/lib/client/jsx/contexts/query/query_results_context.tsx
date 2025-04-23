import React, {useState, createContext, useContext, useCallback} from 'react';

import {QueryResponse} from './query_types';
import {QueryGraphContext} from '../../contexts/query/query_graph_context';
import {QueryColumnContext} from '../../contexts/query/query_column_context';
import {QueryWhereContext} from '../../contexts/query/query_where_context';
import {unpackParams} from '../../utils/query_uri_params';

export const defaultQueryResultsParams = {
  expandMatrices: true,
  flattenQuery: true,
  showDisconnected: false,
  page: 0,
  pageSize: 10,
  data: {} as QueryResponse,
  numRecords: 0,
  queryString: '',
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
    state: {rootModel=''},
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
        setResultsState
      }}
    >
      {props.children}
    </QueryResultsContext.Provider>
  );
};
