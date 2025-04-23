import React, {useEffect, useContext} from 'react';
import {QueryColumnState} from './query_column_context';
import {QueryWhereState} from './query_where_context';
import {QueryColumn} from './query_types';
import {migrateSubclauses, migrateSlices} from '../../utils/query_uri_params';
import {QueryGraphContext} from '../../contexts/query/query_graph_context';
import {QueryColumnContext} from '../../contexts/query/query_column_context';
import {QueryWhereContext} from '../../contexts/query/query_where_context';
import {QueryResultsContext} from '../../contexts/query/query_results_context';
import {packParams, unpackParams} from '../../utils/query_uri_params';

export default function useUriQueryParams({
  search,
  pathname
}: {
  search?: string;
  pathname?: string;
}) {
  if (!search) search = window.location.hash || '';
  if (!pathname) pathname = window.location.pathname;

  const {state: {rootModel=''}} = useContext(QueryGraphContext);
  const {state: columnState} = useContext(QueryColumnContext);
  const {state: whereState} = useContext(QueryWhereContext);
  const {setQueryStateFromString} = useContext(QueryResultsContext);

  // Update the search params to reflect current state
  useEffect(() => {
    const updateSearchParams = async () => {
      let params = await unpackParams(search!.slice(1));
      let paramString = await packParams({
        ...params,
        rootModel,
        ...whereState,
        ...columnState
      });
      if (search!.slice(1) !== paramString) {
        history.pushState({}, '', `${pathname}#${paramString}`);
      }
    }

    updateSearchParams().catch(console.error);
  }, [rootModel, search, whereState, columnState, pathname]);

  // Set current state to reflect query params only on component load
  useEffect(() => {
    const updateState = async () => {
      if (search === '') return;

      const paramString = await packParams({rootModel, ...whereState, ...columnState});
      let serializedState = '#' + paramString;

      if (serializedState === search) return;
      await setQueryStateFromString(search.slice(1));
    }

    updateState().catch(console.error);
  }, []);
}
