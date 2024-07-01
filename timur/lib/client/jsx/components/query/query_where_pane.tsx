import React, {useMemo, useContext, useState, useCallback} from 'react';

import Grid from '@material-ui/core/Grid';
import Button from '@material-ui/core/Button';
import Typography from '@material-ui/core/Typography';
import AddIcon from '@material-ui/icons/Add';

import {makeStyles} from '@material-ui/core/styles';

import {QueryGraphContext} from '../../contexts/query/query_graph_context';
import {QueryWhereContext} from '../../contexts/query/query_where_context';
import {
  QueryFilter,
  QuerySlice,
  EmptyQueryClause
} from '../../contexts/query/query_types';
import QueryClause from './query_clause';
import QueryFilterControl from './query_filter_control';

const useStyles = makeStyles((theme) => ({
  header: {
    marginBottom: 10,
    paddingLeft: '1rem',
    paddingRight: '1rem'
  },
  clauseHeading: {
    paddingLeft: '0.25rem'
  },
  clauseSubheading: {
    color: 'gray',
    fontSize: '0.9rem'
  },
  shimLeft: {
    marginLeft: '-3px'
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

  if (!rootModel) return null;

  return (
    <QueryClause title=''>
      Where:
      {!recordFilters.length && <Typography onClick={addNewRecordFilter} color='gray'>no conditions</Typography>}
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
      {recordFilters.length && <Button onClick={addNewRecordFilter} startIcon={<AddIcon />}>
        Filter
      </Button>
      }
    </QueryClause>
  );
};

export default QueryWherePane;
