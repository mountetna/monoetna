import React from 'react';

const queryString = require('query-string');

import useUriQueryParams from '../use_uri_query_params';

function convertToQueryString(hash: {[key: string]: any}): string {
  return queryString.stringify(
    Object.entries(hash).reduce((acc: any, [key, value]) => {
      acc[key] = JSON.stringify(value);
      return acc;
    }, {})
  );
}

const MockComponent = ({
  setWhereState,
  setQueryColumns,
  setRootModel,
  input
}: {
  setWhereState: () => void;
  setQueryColumns: () => void;
  setRootModel: () => void;
  input: {[key: string]: any};
}) => {
  useUriQueryParams({
    rootModel: 'labor',
    columnState: {
      columns: []
    },
    whereState: {
      recordFilters: [],
      orRecordFilterIndices: []
    },
    setWhereState,
    setQueryColumns,
    setRootModel,
    search: convertToQueryString(input)
  });

  return <>Rendered</>;
};

export default MockComponent;
