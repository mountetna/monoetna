// Generic filter component?
// Model, attribute, operator, operand

import React, {useCallback} from 'react';
import Grid from '@material-ui/core/Grid';
import Button from '@material-ui/core/Button';
import {makeStyles} from '@material-ui/core/styles';
import Tooltip from '@material-ui/core/Tooltip';
import AddIcon from '@material-ui/icons/Add';

import {
  EmptyQuerySubclause,
  QueryClause,
  QuerySubclause
} from '../../contexts/query/query_types';
import {emptyQueryClauseStamp} from '../../selectors/query_selector';
import useQueryClause from './query_use_query_clause';
import {QueryGraph} from '../../utils/query_graph';
import RemoveIcon from './query_remove_icon';
import Selector from './query_selector';
import QueryFilterSubclause from './query_filter_subclause';

const useStyles = makeStyles((theme) => ({
  addSubclauseBtn: {
    paddingLeft: 'calc(0.5rem - 4px)',
    marginBottom: '0.5rem'
  },
  grid: {
    paddingTop: '0.5rem'
  },
  fullWidth: {
    width: '100%'
  }
}));

const QueryFilterClause = ({
  clause,
  modelNames,
  graph,
  isColumnFilter,
  waitTime,
  eager,
  showRemoveIcon,
  canAddSubclause,
  patchClause,
  removeClause
}: {
  clause: QueryClause;
  modelNames: string[];
  graph: QueryGraph;
  isColumnFilter: boolean;
  waitTime?: number;
  eager?: boolean;
  showRemoveIcon: boolean;
  canAddSubclause?: boolean;
  patchClause: (clause: QueryClause) => void;
  removeClause: () => void;
}) => {
  const {modelAttributes} = useQueryClause({
    modelName: clause.modelName,
    graph,
    isColumnFilter
  });
  const classes = useStyles();

  const handleModelSelect = useCallback(
    (modelName: string) => {
      patchClause(emptyQueryClauseStamp(modelName));
    },
    [patchClause]
  );

  const handleUpdateSubclause = useCallback(
    (updatedSubclause: QuerySubclause, index: number) => {
      let updatedSubclauses = [...(clause.subclauses || [])];
      updatedSubclauses[index] = {...updatedSubclause};

      patchClause({
        ...clause,
        subclauses: [...updatedSubclauses]
      });
    },
    [clause, patchClause]
  );

  const handleRemoveSubclause = useCallback(
    (index: number) => {
      let updatedSubclauses = [...(clause.subclauses || [])];
      updatedSubclauses.splice(index, 1);
      if (0 === updatedSubclauses.length) {
        updatedSubclauses.push({...EmptyQuerySubclause});
      }
      patchClause({
        ...clause,
        subclauses: [...updatedSubclauses]
      });
    },
    [clause, patchClause]
  );

  const handleAddSubclause = useCallback(() => {
    let updatedSubclauses = [...(clause.subclauses || [])];
    updatedSubclauses.push({...EmptyQuerySubclause});
    patchClause({
      ...clause,
      subclauses: [...updatedSubclauses]
    });
  }, [clause, patchClause]);

  return (
    <Grid container alignItems='center'>
      <Grid item xs={3}>
        <Selector
          canEdit={true}
          name={clause.modelName}
          onSelect={handleModelSelect}
          choiceSet={modelNames}
          label='model'
        />
      </Grid>
      <Grid
        item
        xs={showRemoveIcon ? 8 : 9}
        container
        direction='column'
        alignItems='center'
        className={classes.fullWidth}
      >
        <Grid item className={classes.fullWidth}>
          {clause.subclauses?.map(
            (subclause: QuerySubclause, index: number) => {
              return (
                <QueryFilterSubclause
                  subclause={subclause}
                  subclauseIndex={index}
                  graph={graph}
                  waitTime={waitTime}
                  eager={eager}
                  modelName={clause.modelName}
                  modelAttributes={modelAttributes}
                  patchSubclause={(updatedSubclause: QuerySubclause) =>
                    handleUpdateSubclause(updatedSubclause, index)
                  }
                  removeSubclause={() => handleRemoveSubclause(index)}
                  isColumnFilter={isColumnFilter}
                  showRemoveIcon={!isColumnFilter}
                />
              );
            }
          )}
        </Grid>
        <Grid item className={classes.fullWidth}>
          {canAddSubclause ? (
            <Tooltip title='Add subclause' aria-label='Add subclause'>
              <Button
                className={classes.addSubclauseBtn}
                startIcon={<AddIcon />}
                onClick={handleAddSubclause}
              >
                Subclause
              </Button>
            </Tooltip>
          ) : null}
        </Grid>
      </Grid>
      {showRemoveIcon ? (
        <Grid item xs={1} alignItems='center'>
          <RemoveIcon
            showRemoveIcon={showRemoveIcon}
            onClick={removeClause}
            label='clause'
          />
        </Grid>
      ) : null}
    </Grid>
  );
};
export default QueryFilterClause;
