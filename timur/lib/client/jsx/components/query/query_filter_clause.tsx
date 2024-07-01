// Generic filter component?
// Model, attribute, operator, operand

import React, {useCallback} from 'react';
import Grid from '@material-ui/core/Grid';
import Button from '@material-ui/core/Button';
import {makeStyles} from '@material-ui/core/styles';
import Tooltip from '@material-ui/core/Tooltip';
import Typography from '@material-ui/core/Typography';
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
import QueryNumber from './query_number';
import Selector from './query_selector';
import QueryFilterSubclause from './query_filter_subclause';

const useStyles = makeStyles((theme) => ({
  filter_clause: {
    paddingLeft: '25px'
  },
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
  canAddSubclause = false,
  hasModelChildren,
  patchClause,
  index,
  removeClause
}: {
  clause: QueryClause;
  modelNames: string[];
  graph: QueryGraph;
  isColumnFilter: boolean;
  waitTime?: number;
  eager?: boolean;
  index?: number;
  showRemoveIcon: boolean;
  canAddSubclause: boolean;
  hasModelChildren: boolean;
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
     <Grid item container alignItems='center'>
        <RemoveIcon
          showRemoveIcon={showRemoveIcon}
          onClick={removeClause}
          label='clause'
        />
        <QueryNumber number={index} level={1}/>
        {
          hasModelChildren && <Selector
            canEdit={true}
            name={clause.any ? 'Any' : 'Every'}
            onSelect={(val: string) =>
              handleClauseAnySelect(val, clause, index)
            }
            choiceSet={['Any', 'Every']}
            label='model'
          />
        }
        <MapSelector
          canEdit={true}
          name={clause.modelName}
          onSelect={handleModelSelect}
          choiceSet={modelNames}
          label='model'
        />
        <Grid container className={classes.filter_clause} direction='column' alignItems='flex-start'>
          {clause.subclauses?.map(
            (subclause: QuerySubclause, index: number) => {
              return (
                <QueryFilterSubclause
                  key={index}
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
                  showRemoveIcon={canAddSubclause}
                />
              );
            }
          )}
        </Grid>
        <Grid item className={classes.fullWidth}>
          {canAddSubclause ? (
            <Tooltip title='Add subclause' aria-label='Add subclause'>
              <Button
                variant='text'
                className={classes.addSubclauseBtn}
                startIcon={<AddIcon />}
                onClick={handleAddSubclause}
              >
                <QueryNumber number={clause.subclauses.length} level={2}/>
                Add subclause
              </Button>
            </Tooltip>
          ) : null}
      </Grid>
    </Grid>
  );
};
export default QueryFilterClause;
