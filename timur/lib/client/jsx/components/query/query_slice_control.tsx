import React, {useCallback} from 'react';
import Grid from '@material-ui/core/Grid';

import {QueryClause, QuerySlice} from '../../contexts/query/query_types';
import {QueryGraph} from '../../utils/query_graph';
import QueryFilterClause from './query_filter_clause';
import RemoveIcon from './query_remove_icon';
import QueryNumber from './query_number';

const QuerySliceControl = ({
  slice,
  modelNames,
  graph,
  patchSlice,
  removeSlice
}: {
  slice: QuerySlice;
  modelNames: string[];
  graph: QueryGraph;
  patchSlice: (slice: QuerySlice) => void;
  removeSlice: () => void;
}) => {
  const handlePatchClause = useCallback(
    (clause: QueryClause) => {
      patchSlice({
        ...slice,
        clause
      });
    },
    [patchSlice, slice]
  );

  return <Grid item container alignItems='center'>
    <QueryFilterClause
      index={sliceIndex}
      clause={slice.clause}
      graph={graph}
      modelNames={modelNames}
      isColumnFilter={true}
      patchClause={handlePatchClause}
      removeClause={removeSlice}
      showRemoveIcon={true}
      canAddSubclause={false}
    />
  </Grid>;
};
export default QuerySliceControl;
