import React, {useCallback} from 'react';
import Grid from '@material-ui/core/Grid';

import {QueryClause, QuerySlice} from '../../contexts/query/query_types';
import {QueryGraph} from '../../utils/query_graph';
import QueryFilterClause from './query_filter_clause';
import RemoveIcon from './query_remove_icon';

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

  return (
    <>
      <Grid item container xs={11}>
        <QueryFilterClause
          clause={slice.clause}
          graph={graph}
          modelNames={modelNames}
          isColumnFilter={true}
          patchClause={handlePatchClause}
          removeClause={() => {}}
          showRemoveIcon={false}
          canAddSubclause={false}
        />
      </Grid>
      <Grid item xs={1} container justify='flex-end'>
        <RemoveIcon showRemoveIcon={true} onClick={removeSlice} label='slice' />
      </Grid>
    </>
  );
};
export default QuerySliceControl;
