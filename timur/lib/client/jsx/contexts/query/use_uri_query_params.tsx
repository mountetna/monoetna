import React, {useEffect} from 'react';
import {QueryColumnState} from './query_column_context';
import {QueryWhereState} from './query_where_context';
import {QueryColumn} from './query_types';
import {migrateSubclauses, migrateSlices} from '../../utils/query_uri_params';

const btou = (str: string) => btoa(str).replace(/\+/g, '-').replace(/\//g, '_').replace(/=/g, '');
const utob = (str: string) => {
  str = str.replace(/-/g, '+').replace(/_/g, '/');
  while (str.length % 4) str += '=';
  return atob(str);
}

const unpackParams = async (search: string) => {
  let searchParams = new URLSearchParams(search);

  if (searchParams.has('q')) {
    const stream = new Blob(
      [Uint8Array.from(
        utob(searchParams.get('q') as string),
        c => c.charCodeAt(0)
      )]
    ).stream();
    let decompressedStream = stream.pipeThrough(
      new DecompressionStream("gzip")
    );
    const blob = await new Response(decompressedStream as BodyInit).blob();
    const text = await blob.text();
    return JSON.parse(text);
  }

  return {
    rootModel: searchParams.get('rootModel') || '',
    recordFilters: migrateSubclauses(
      JSON.parse(searchParams.get('recordFilters') || '[]')
    ),
    orRecordFilterIndices: JSON.parse(
      searchParams.get('orRecordFilterIndices') || '[]'
    ),
    columns: migrateSlices(JSON.parse(searchParams.get('columns') || '[]'))
  }
};

const packParams = async (params: any) => {
  const stream = new Blob(
    [JSON.stringify(params)],
    { type: 'applicaton/json' }
  ).stream();

  const gzipStream = stream.pipeThrough(new CompressionStream("gzip"));

  const compressedResponse = await new Response(gzipStream as BodyInit);
  const blob = await compressedResponse.blob();
  const buffer = await blob.arrayBuffer();
  const compressedBase64 = btou(
    String.fromCharCode(
      ...new Uint8Array(buffer)
    )
  );

  const searchParams = new URLSearchParams();

  searchParams.set('q', compressedBase64)

  return searchParams.toString();
}

function serializeState(state: {[key: string]: any}, isJson: boolean = true) {
  let params = new URLSearchParams();
  Object.entries(state).forEach(([key, value]) => {
    params.set(key, (typeof value  === 'string') ? value : JSON.stringify(value));
  })
  return params.toString();
}

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
  if (!search) search = window.location.hash || '';
  if (!pathname) pathname = window.location.pathname;

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
    const setCurrentState = async () => {
      if (search === '') return;

      const paramString = await packParams({rootModel, ...whereState, ...columnState});
      let serializedState = '#' + paramString;

      if (serializedState === search) return;

      let params = await unpackParams(search.slice(1));

      setQueryColumns(params.columns);
      setRootModel(params.rootModel);
      setWhereState({
        recordFilters: params.recordFilters,
        orRecordFilterIndices: params.orRecordFilterIndices
      });
    }

    setCurrentState().catch(console.error);

  }, []);
}
