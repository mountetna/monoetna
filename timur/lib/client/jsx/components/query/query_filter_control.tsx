import React, {useCallback, useMemo, useState, useContext} from 'react';
import Grid from '@material-ui/core/Grid';
import IconButton from '@material-ui/core/IconButton';
import {makeStyles} from '@material-ui/core/styles';
import Tooltip from '@material-ui/core/Tooltip';
import Checkbox from '@material-ui/core/Checkbox';
import Typography from '@material-ui/core/Typography';
import AddIcon from '@material-ui/icons/Add';
import FileCopyIcon from '@material-ui/icons/FileCopy';

import {QueryClause, QueryFilter} from '../../contexts/query/query_types';
import {emptyQueryClauseStamp} from '../../selectors/query_selector';
import {QueryGraph} from '../../utils/query_graph';
import {QueryGraphContext} from '../../contexts/query/query_graph_context';
import QueryFilterClause from './query_filter_clause';
import QueryNumber from './query_number';
import QueryAnyEverySelectorList from './query_any_every_selector_list';
import RemoveIcon from './query_remove_icon';
import QueryModelSelector from './query_model_selector';

const useStyles = makeStyles((theme) => ({
  paper: {
    padding: '0.5rem 0.5rem 0 0.5rem',
    marginBottom: '0.5rem',
    minHeight: '48px'
  },
  paddingLeft: {
    paddingLeft: 'calc(0.5rem - 4px)'
  },
  and: {
    padding: '10px 0px',
    cursor: 'pointer'
  },
  grid: {
    paddingLeft: '25px'
  }
}));

const QueryFilterControl = ({
  filter,
  modelNames,
  or,
  setOr,
  filterIndex,
  patchRecordFilter,
  patchFilter,
  removeFilter,
  copyFilter
}: {
  filter: QueryFilter;
  patchRecordFilter: (index: number, updatedFilter: QueryFilter) => void;
  or: boolean;
  setOr: () => void;
  filterIndex: number;
  modelNames: string[];
  patchFilter: (filter: QueryFilter) => void;
  removeFilter: () => void;
  copyFilter: () => void;
}) => {
  const classes = useStyles();

  const { state: {graph} } = useContext(QueryGraphContext);

  const handleModelSelect = useCallback(
    (modelName: string) => {
      patchFilter({
        modelName,
        anyMap: {},
        clauses: [emptyQueryClauseStamp(modelName)]
      });
    },
    [patchFilter]
  );

  const handlePatchClause = useCallback(
    (clause: QueryClause, index: number) => {
      let updatedClauses = [...filter.clauses];
      updatedClauses[index] = clause;
      patchFilter({
        ...filter,
        clauses: updatedClauses
      });
    },
    [patchFilter, filter]
  );

  const handleRemoveClause = useCallback(
    (index: number) => {
      let updatedClauses = [...filter.clauses];
      updatedClauses.splice(index, 1);
      patchFilter({
        ...filter,
        clauses: updatedClauses
      });
    },
    [patchFilter, filter]
  );

  const handleAddClause = useCallback(() => {
    patchFilter({
      ...filter,
      clauses: [...filter.clauses].concat([
        emptyQueryClauseStamp(filter.modelName)
      ])
    });
  }, [patchFilter, filter]);

  const handleClauseAnySelect = useCallback(
    (val: string, clause: QueryClause, index: number) => {
      handlePatchClause(
        {
          ...clause,
          any: val === 'Any'
        },
        index
      );
    },
    [handlePatchClause]
  );

  const filterModelsMap = graph.subgraphMap(filter.modelName);

  const [ removeHint, setRemoveHint ] = useState(false);
  const [ showControls, setShowControls ] = useState(false);

  return <Grid
    style={{ textDecoration: removeHint ? 'line-through' : 'none' }}
    container
    alignItems='center'
    justifyContent='flex-start'
    className='query-where-selector'
    onMouseEnter={ () => setShowControls(true) }
    onMouseLeave={ () => setShowControls(false) }
  >
    { filterIndex > 0 && <Grid className={classes.and} container>
    <Typography style={{ color:'purple'}} onClick={setOr}>{ or ? 'or' : 'and' }</Typography>
    </Grid>
    }
    <QueryNumber setRemoveHint={ setRemoveHint } onClick={ removeFilter } number={filterIndex} level={0}/>
    <QueryAnyEverySelectorList
      filter={filter}
      index={filterIndex}
      patchRecordFilter={patchRecordFilter}
    />
    <QueryModelSelector
      modelName={filter.modelName}
      setModel={handleModelSelect}
      modelNames={modelNames}
    />
    { showControls && <>
        <Tooltip title='Add clause' aria-label='Add clause'>
          <IconButton size='small' onClick={handleAddClause} color='primary'>
            <AddIcon fontSize='small'/>
          </IconButton>
        </Tooltip>
        <Tooltip title='Copy filter' aria-label='Copy filter'>
          <IconButton size='small' onClick={copyFilter} color='primary'>
            <FileCopyIcon fontSize='small'/>
          </IconButton>
        </Tooltip>
      </>
      }
    <Grid item container direction='column'>
      <Grid container direction='column' className={classes.grid} alignItems='flex-start'>
        {filter.clauses.map((clause: QueryClause, index: number) => {
          return (
              <QueryFilterClause
                key={index}
                index={index}
                clause={clause}
                modelNames={Object.keys(filterModelsMap)}
                hasModelChildren={!!filterModelsMap[clause.modelName]}
                patchClause={(updatedClause: QueryClause) =>
                  handlePatchClause(updatedClause, index)
                }
                selectClause={ (val: string) =>
                  handleClauseAnySelect(val, clause, index)
                }
                removeClause={() => handleRemoveClause(index)}
                showRemoveIcon={
                  !(0 === index && 1 === filter.clauses.length)
                }
                canAddSubclause={true}
              />
          );
        })}
      </Grid>
    </Grid>
  </Grid>;
};
export default QueryFilterControl;
