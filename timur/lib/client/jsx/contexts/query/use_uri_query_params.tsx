import React, {useEffect} from 'react';
import {QueryColumnState} from './query_column_context';
import {QueryWhereState} from './query_where_context';
import {QueryColumn} from './query_types';
import {migrateSubclauses, migrateSlices} from '../../utils/query_uri_params';

export default function useUriQueryParams({
  columnState,
  rootModel,
  whereState,
  setQueryColumns,
  setRootModel,
  setWhereState,
  search,
  pathname
}: {
  columnState: QueryColumnState;
  rootModel: string;
  whereState: QueryWhereState;
  setQueryColumns: (columns: QueryColumn[]) => void;
  setRootModel: (modelName: string) => void;
  setWhereState: (whereState: QueryWhereState) => void;
  search?: string;
  pathname?: string;
}) {
  if (!search) search = window.location.hash;
  if (!pathname) pathname = window.location.pathname;

  function serializeState(state: {[key: string]: any}, isJson: boolean = true) {
    return Object.entries(state)
      .map(([key, value]) => {
        return `${key}=${encodeURIComponent(
          isJson ? JSON.stringify(value) : value
        )}`;
      })
      .join('&');
  }

  // Update the search params to reflect current state
  useEffect(() => {
    let searchParams = new URLSearchParams(search.slice(1));
    searchParams.set('rootModel', rootModel);
    Object.entries(whereState).forEach(([key, value]: [string, any]) => {
      searchParams.set(key, JSON.stringify(value));
    });
    Object.entries(columnState).forEach(([key, value]: [string, any]) => {
      searchParams.set(key, JSON.stringify(value));
    });

    if (search.slice(1) !== searchParams.toString()) {
      history.pushState({}, '', `${pathname}#${searchParams.toString()}`);
    }
  }, [rootModel, search, whereState, columnState, pathname]);

  // Set current state to reflect query params only on component load
  useEffect(() => {
    if (search === '') return;

    let searchParams = new URLSearchParams(search.slice(1));

    let serializedState =
      '#' +
      serializeState({rootModel}, false) +
      serializeState(whereState) +
      serializeState(columnState);

    if (serializedState === search) return;

    setQueryColumns(
      migrateSlices(JSON.parse(searchParams.get('columns') || '[]'))
    );
    setRootModel(searchParams.get('rootModel') || '');
    setWhereState({
      recordFilters: migrateSubclauses(
        JSON.parse(searchParams.get('recordFilters') || '[]')
      ),
      orRecordFilterIndices: JSON.parse(
        searchParams.get('orRecordFilterIndices') || '[]'
      )
    });
  }, []);
}
