// Generic filter component?
// Model, attribute, operator, operand

import React, {useState,useCallback} from 'react';
import Grid from '@material-ui/core/Grid';
import Button from '@material-ui/core/Button';
import {makeStyles} from '@material-ui/core/styles';
import Tooltip from '@material-ui/core/Tooltip';
import IconButton from '@material-ui/core/IconButton';
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
import QueryModelSelector from './query_model_selector';
import QueryFilterSubclause from './query_filter_subclause';

const useStyles = makeStyles((theme) => ({
  filter_clause: {
    paddingLeft: '25px'
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
  selectClause,
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
  selectClause: (val: string) => void;
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

  const [ removeHint, setRemoveHint ] = useState(false);
  const [ showControls, setShowControls ] = useState(false);

  return (
    <Grid
      style={{ textDecoration: removeHint ? 'line-through' : 'none' }}
      onMouseEnter={ () => setShowControls(true) }
      onMouseLeave={ () => setShowControls(false) }
      item container alignItems='center'>
      <QueryNumber
        setRemoveHint={ showRemoveIcon ? setRemoveHint : null }
        onClick={ showRemoveIcon ? removeClause : null}
        number={index} level={1}/>
        {
          hasModelChildren && <Selector
            canEdit={true}
            name={clause.any ? 'Any' : 'Every'}
            onSelect={selectClause}
            choiceSet={['Any', 'Every']}
            label='model'
          />
        }
        <QueryModelSelector
          modelName={clause.modelName}
          setModel={handleModelSelect}
          options={modelNames}
        />
          {canAddSubclause && showControls && (
            <Tooltip title='Add subclause' aria-label='Add subclause'>
              <IconButton size='small' onClick={handleAddSubclause} color='primary'>
                <AddIcon fontSize='small'/>
              </IconButton>
            </Tooltip>
          )}
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
    </Grid>
  );
};
export default QueryFilterClause;
