import React, {useCallback, useMemo} from 'react';
import Grid from '@material-ui/core/Grid';
import Button from '@material-ui/core/Button';
import {makeStyles} from '@material-ui/core/styles';
import Tooltip from '@material-ui/core/Tooltip';
import Checkbox from '@material-ui/core/Checkbox';
import Typography from '@material-ui/core/Typography';
import AddIcon from '@material-ui/icons/Add';

import {QueryClause, QueryFilter} from '../../contexts/query/query_types';
import {emptyQueryClauseStamp} from '../../selectors/query_selector';
import {QueryGraph} from '../../utils/query_graph';
import QueryFilterClause from './query_filter_clause';
import QueryNumber from './query_number';
import QueryAnyEverySelectorList from './query_any_every_selector_list';
import RemoveIcon from './query_remove_icon';
import CopyIcon from './query_copy_icon';
import Selector from './query_selector';
import MapSelector from './map_selector';

const useStyles = makeStyles((theme) => ({
  paper: {
    padding: '0.5rem 0.5rem 0 0.5rem',
    marginBottom: '0.5rem',
    minHeight: '48px'
  },
  paddingLeft: {
    paddingLeft: 'calc(0.5rem - 4px)'
  },
  grid: {
    paddingTop: '0.5rem',
    paddingLeft: '25px'
  }
}));

const QueryFilterControl = ({
  filter,
  modelNames,
  graph,
  or,
  setOr,
  filterIndex,
  patchRecordFilter,
  patchFilter,
  removeFilter,
  copyFilter
}: {
  filter: QueryFilter;
  patchRecordFilter: QueryFilter;
  modelNames: string[];
  graph: QueryGraph;
  patchFilter: (filter: QueryFilter) => void;
  removeFilter: () => void;
  copyFilter: () => void;
}) => {
  const classes = useStyles();

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

  const modelChildren = useMemo(() => {
    if (!filter.modelName || filter.modelName === '') return {};

    return graph.childrenMap(filter.modelName);
  }, [filter.modelName, graph]);

  return <Grid
    container
    alignItems='center'
    justifyContent='flex-start'
    className='query-where-selector'
  >
    <CopyIcon canEdit={true} onClick={copyFilter} label='filter' />
    <RemoveIcon
      showRemoveIcon={true}
      onClick={removeFilter}
      label='filter'
    />
    <QueryNumber number={filterIndex} level={0}/>
    <Checkbox
      checked={or}
      color='primary'
      onChange={setOr}
      inputProps={{'aria-label': 'secondary checkbox'}}
    />
    <QueryAnyEverySelectorList
      filter={filter}
      index={filterIndex}
      patchRecordFilter={patchRecordFilter}
    />
    <MapSelector
      canEdit={true}
      modelName={filter.modelName}
      setModel={handleModelSelect}
      options={modelNames}
    />
    <Grid item container direction='column'>
      <Grid container direction='column' className={classes.grid} alignItems='flex-start'>
        {filter.clauses.map((clause: QueryClause, index: number) => {
          return (
              <QueryFilterClause
                key={index}
                index={index}
                clause={clause}
                graph={graph}
                modelNames={Object.keys(modelChildren)}
                hasModelChildren={!!modelChildren[clause.modelName]}
                isColumnFilter={false}
                patchClause={(updatedClause: QueryClause) =>
                  handlePatchClause(updatedClause, index)
                }
                removeClause={() => handleRemoveClause(index)}
                showRemoveIcon={
                  !(0 === index && 1 === filter.clauses.length)
                }
                canAddSubclause={true}
              />
          );
        })}
        <Tooltip title='Add and clause' aria-label='Add and clause'>
          <Button
            variant='text'
            className={classes.paddingLeft}
            startIcon={<AddIcon />}
            onClick={handleAddClause}
          >
            <QueryNumber number={filter.clauses.length} level={1}/>
            Add clause
          </Button>
        </Tooltip>
      </Grid>
    </Grid>
  </Grid>;
};
export default QueryFilterControl;
