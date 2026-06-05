import React, {useState, createContext, useCallback} from 'react';

import {QueryFilter} from './query_types';

export const defaultQueryWhereParams = {
  recordFilters: [] as QueryFilter[],
  orRecordFilterIndices: [] as number[],
  globalOr: false
};

const defaultQueryWhereState = {
  ...defaultQueryWhereParams
};

export type QueryWhereState = Readonly<typeof defaultQueryWhereState>;

export const defaultQueryWhereContext = {
  state: defaultQueryWhereState as QueryWhereState,
  addRecordFilter: (recordFilter: QueryFilter) => {},
  removeRecordFilter: (index: number) => {},
  removeAllRecordFilters: () => {},
  patchRecordFilter: (index: number, recordFilter: QueryFilter) => {},
  setOrRecordFilterIndices: (indices: number[]) => {},
  setGlobalOr: (or: boolean) => {},
  setWhereState: (newState: QueryWhereState) => {},
  resetWhereState: () => {}
};

export type QueryWhereContextData = typeof defaultQueryWhereContext;

export const QueryWhereContext = createContext(defaultQueryWhereContext);
export type QueryWhereContext = typeof QueryWhereContext;
export type ProviderProps = {params?: {}; children: any};

export const QueryWhereProvider = (
  props: ProviderProps & Partial<QueryWhereContextData>
) => {
  const [state, setState] = useState(
    props.state || defaultQueryWhereContext.state
  );

  const addRecordFilter = useCallback(
    (recordFilter: QueryFilter) => {
      setState({
        ...state,
        recordFilters: [...(state.recordFilters || [])].concat([recordFilter])
      });
    },
    [state]
  );

  const removeRecordFilter = useCallback(
    (filterIndex: number) => {
      let updatedRecordFilters = [...state.recordFilters];
      let updatedOrRecordFilterIndices = state.orRecordFilterIndices
        .filter((index) => index !== filterIndex)
        .map((index) => {
          if (index >= filterIndex) return index - 1;
          return index;
        });
      updatedRecordFilters.splice(filterIndex, 1);
      setState({
        ...state,
        recordFilters: updatedRecordFilters,
        orRecordFilterIndices: updatedOrRecordFilterIndices
      });
    },
    [state]
  );

  const removeAllRecordFilters = useCallback(
    () => {
      setState({
        ...state,
        recordFilters: [],
        orRecordFilterIndices: [],
        globalOr: false
      });
    },
    [state]
  );

  const patchRecordFilter = useCallback(
    (index: number, recordFilter: QueryFilter) => {
      let updatedRecordFilters = [...state.recordFilters];
      updatedRecordFilters[index] = recordFilter;
      setState({
        ...state,
        recordFilters: updatedRecordFilters
      });
    },
    [state]
  );

  const setOrRecordFilterIndices = useCallback(
    (orRecordFilterIndices: number[]) => {
      setState({
        ...state,
        orRecordFilterIndices
      });
    },
    [state]
  );

  const setGlobalOr = useCallback(
    (globalOr: boolean) => {
      setState({
        ...state,
        globalOr
      });
    },
    [state]
  );

  const setWhereState = useCallback((newState: QueryWhereState) => {
    setState({
      ...newState
    });
  }, []);

  const resetWhereState = useCallback(() => {
    setState({
      ...defaultQueryWhereParams
    });
  }, []);

  return (
    <QueryWhereContext.Provider
      value={{
        state,
        addRecordFilter,
        removeRecordFilter,
        removeAllRecordFilters,
        patchRecordFilter,
        setOrRecordFilterIndices,
        setWhereState,
        setGlobalOr,
        resetWhereState
      }}
    >
      {props.children}
    </QueryWhereContext.Provider>
  );
};
