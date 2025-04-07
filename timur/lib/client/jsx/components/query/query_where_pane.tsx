import React, {useMemo, useContext, useState, useCallback} from 'react';

import Grid from '@material-ui/core/Grid';
import IconButton from '@material-ui/core/IconButton';
import Tooltip from '@material-ui/core/Tooltip';
import Typography from '@material-ui/core/Typography';
import AddIcon from '@material-ui/icons/Add';
import RestartAltIcon from '@material-ui/icons/RestartAlt';

import {makeStyles} from '@material-ui/core/styles';

import {QueryGraphContext} from '../../contexts/query/query_graph_context';
import {QueryWhereContext} from '../../contexts/query/query_where_context';
import {
  QueryFilter,
  QuerySlice,
  EmptyQueryClause
} from '../../contexts/query/query_types';
import QueryClause from './query_clause';
import QueryClauseSummaryControls from './query_clause_summary_controls';
import QueryChevron from './query_chevron';
import QueryFilterControl from './query_filter_control';

const useStyles = makeStyles((theme) => ({
  folded: {
    fontStyle: 'italic',
    paddingLeft: '10px',
    cursor: 'pointer'
  },
  conditions: {
    width: '100%',
    paddingLeft: '30px'
  },
  title: {
    fontWeight: 'bold'
  },
  empty: {
    paddingLeft: '5px'
  }
}));

const QueryWherePane = () => {
  // Use an update counter to get the child components
  //  (i.e. the QueryFilterControls) to remount whenever
  //  the recordFilters list has one removed.
  // If not, the component rendered state gets confused
  //  because of non-unique keys.
  const [updateCounter, setUpdateCounter] = useState(0);
  const {
    state: {graph, rootModel}
  } = useContext(QueryGraphContext);
  const {
    state: {orRecordFilterIndices, recordFilters},
    addRecordFilter,
    removeRecordFilter,
    removeAllRecordFilters,
    patchRecordFilter,
    setOrRecordFilterIndices
  } = useContext(QueryWhereContext);
  const classes = useStyles();

  const addNewRecordFilter = useCallback(() => {
    addRecordFilter({
      modelName: '',
      clauses: [
        {
          ...EmptyQueryClause
        }
      ],
      anyMap: {}
    });
  }, [addRecordFilter]);

  const handlePatchFilter = useCallback(
    (
      index: number,
      updatedFilter: QueryFilter,
      originalFilter: QueryFilter
    ) => {
      if (
        updatedFilter.modelName !== originalFilter.modelName &&
        rootModel != null
      ) {
        let selectableModels = graph.sliceableModelNamesInPath(
          rootModel,
          updatedFilter.modelName
        );

        updatedFilter.anyMap = selectableModels.reduce(
          (acc: {[key: string]: boolean}, modelName: string) => {
            acc[modelName] = true;
            return acc;
          },
          {}
        );
      }
      patchRecordFilter(index, updatedFilter);
    },
    [patchRecordFilter, graph, rootModel]
  );

  function handleRemoveFilter(index: number) {
    removeRecordFilter(index);
    setUpdateCounter(updateCounter + 1);
  }

  const handleChangeOrFilters = useCallback(
    (index: number) => {
      let copy = [...orRecordFilterIndices];

      if (copy.includes(index)) copy.splice(copy.indexOf(index), 1);
      else copy.push(index);

      setOrRecordFilterIndices(copy);
    },
    [orRecordFilterIndices, setOrRecordFilterIndices]
  );

  const handleCopyFilter = useCallback(
    (filter: QueryFilter) => {
      addRecordFilter({...filter});
    },
    [addRecordFilter]
  );

  const modelNames = useMemo(
    () => [...new Set(graph.allPaths(rootModel).flat())].sort(),
    [graph, rootModel]
  );

  const [ fold, setFold ] = useState(true);

  if (!rootModel) return null;

  return (
    <QueryClause title=''>
      <Grid container alignItems='center'>
        <QueryChevron fold={fold} setFold={setFold}/>
        <Typography className={classes.title}>Where:</Typography>
        <QueryClauseSummaryControls
          fold={fold}
          setFold={setFold}
          addHandler={addNewRecordFilter}
          removeHandler={removeAllRecordFilters}
          itemName='condition'
          numItems={recordFilters.length}/>
      </Grid>
      {
        !fold && <Grid container direction='column' className={classes.conditions}>
          {!recordFilters.length && <Typography className={classes.empty} onClick={addNewRecordFilter} style={{ color:'gray' }}>no conditions</Typography> }
          {recordFilters.map((filter: QueryFilter, index: number) => (
            <QueryFilterControl
              key={`${index}-${updateCounter}`}
              or={orRecordFilterIndices.includes(index)}
              setOr={(e, checked) => handleChangeOrFilters(index)}
              filterIndex={index}
              filter={filter}
              patchRecordFilter={patchRecordFilter}
              modelNames={modelNames}
              graph={graph}
              patchFilter={(updatedFilter: QueryFilter | QuerySlice) =>
                handlePatchFilter(
                  index,
                  updatedFilter as QueryFilter,
                  filter
                )
              }
              removeFilter={() => handleRemoveFilter(index)}
              copyFilter={() => handleCopyFilter(filter)}
            />
          ))}
        </Grid>
      }
    </QueryClause>
  );
};

export default QueryWherePane;
