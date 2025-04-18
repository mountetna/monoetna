import React, {useCallback, useState, useContext} from 'react';
import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';

import {QueryClause, QuerySubclause, QuerySlice} from '../../contexts/query/query_types';
import {QueryGraph} from '../../utils/query_graph';
import {QueryGraphContext} from '../../contexts/query/query_graph_context';
import RemoveIcon from './query_remove_icon';
import QueryNumber from './query_number';
import QueryModelAttributeSelector from './query_model_attribute_selector';
import QueryFilterSubclause from './query_filter_subclause';
import { makeStyles } from '@material-ui/core/styles';
import {emptyQueryClauseStamp} from '../../selectors/query_selector';

const useStyles = makeStyles((theme) => ({
  subclause: {
    flex: '1',
    paddingLeft: '5px'
  }
}));

const QuerySliceControl = ({
  slice,
  sliceIndex,
  modelNames,
  patchSlice,
  removeSlice
}: {
  slice: QuerySlice;
  sliceIndex: number;
  modelNames: string[];
  patchSlice: (slice: QuerySlice) => void;
  removeSlice: () => void;
}) => {
  const classes = useStyles();

  const { state: {graph} } = useContext(QueryGraphContext);

  const clause = slice.clause;

  const handlePatchClause = useCallback(
    (clause: QueryClause) => {
      patchSlice({
        ...slice,
        clause
      });
    },
    [patchSlice, slice]
  );
  
  const handleModelSelect = useCallback(
    (modelName: string) => {
      handlePatchClause(emptyQueryClauseStamp(modelName));
    },
    [handlePatchClause]
  );

  const handleUpdateSubclause = useCallback(
    (updatedSubclause: QuerySubclause, index: number) => {
      let updatedSubclauses = [...(clause.subclauses || [])];
      updatedSubclauses[index] = {...updatedSubclause};

      handlePatchClause({
        ...clause,
        subclauses: [...updatedSubclauses]
      });
    },
    [clause, handlePatchClause]
  );

  const handleAttributeSelect = useCallback(
    (modelName:string, attributeName: string) => {
      handlePatchClause({
        ...clause,
        modelName,
        subclauses: [
          {
            ...clause.subclauses?.[0],
            attributeName,
            attributeType: '',
            operator: '',
            operand: ''
          }
        ]
      });
    },
    [clause, handlePatchClause]
  );

  const [ removeHint, setRemoveHint ] = useState(false);

  return <Grid item container alignItems='center'>
    <Grid
      style={{ textDecoration: removeHint ? 'line-through' : 'none' }}
      item container alignItems='center'>
      <QueryNumber
        setRemoveHint={ setRemoveHint }
        onClick={ removeSlice }
        number={sliceIndex} level={1}/>
      <QueryModelAttributeSelector
        modelName={clause.modelName}
        attributeName={clause.subclauses?.[0]?.attributeName}
        setModel={handleModelSelect}
        setAttribute={handleAttributeSelect}
        modelNames={modelNames}
      />
      <Grid item container className={classes.subclause}>
        <QueryFilterSubclause
          key={sliceIndex}
          subclause={clause.subclauses?.[0] as QuerySubclause}
          subclauseIndex={0}
          modelName={clause.modelName}
          patchSubclause={(updatedSubclause: QuerySubclause) =>
            handleUpdateSubclause(updatedSubclause, 0)
          }
          removeSubclause={() => {} }
          isColumnFilter={true}
          showRemoveIcon={false}
        />
      </Grid>
    </Grid>
  </Grid>;
};
export default QuerySliceControl;
