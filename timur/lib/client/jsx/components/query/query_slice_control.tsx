import React, {useCallback} from 'react';
import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';

import {QueryClause, QuerySlice} from '../../contexts/query/query_types';
import {QueryGraph} from '../../utils/query_graph';
import QueryFilterClause from './query_filter_clause';
import RemoveIcon from './query_remove_icon';
import QueryNumber from './query_number';
import QueryModelAttributeSelector from './query_model_attribute_selector';
import QueryFilterSubclause from './query_filter_subclause';
import { makeStyles } from '@mui/styles';

const useStyles = makeStyles((theme) => ({
  subclause: {
    flex: '1',
    paddingLeft: '5px'
  }
}));

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
  const classes = useStyles();

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
            ...clause.subclauses[0],
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


    const x = <QueryFilterClause
      index={sliceIndex}
      clause={slice.clause}
      graph={graph}
      modelNames={modelNames}
      isColumnFilter={true}
      patchClause={handlePatchClause}
      removeClause={removeSlice}
      showRemoveIcon={true}
      canAddSubclause={false}
    />;

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
        attributeName={clause.subclauses[0].attributeName}
        setModel={handleModelSelect}
        setAttribute={handleAttributeSelect}
        modelNames={modelNames}
      />
      <Grid item container className={classes.subclause}>
        <QueryFilterSubclause
          key={sliceIndex}
          subclause={clause.subclauses[0]}
          subclauseIndex={sliceIndex}
          graph={graph}
          modelName={clause.modelName}
          patchSubclause={(updatedSubclause: QuerySubclause) =>
            handleUpdateSubclause(updatedSubclause, sliceIndex)
          }
          modelAttributes={[]}
          removeSubclause={() => handleRemoveSubclause(sliceIndex)}
          isColumnFilter={true}
          showRemoveIcon={false}
        />
      </Grid>
    </Grid>
  </Grid>;
};
export default QuerySliceControl;
