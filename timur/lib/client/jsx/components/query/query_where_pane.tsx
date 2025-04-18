import React, {useMemo, useContext, useState, useCallback} from 'react';

import Checkbox from '@material-ui/core/Checkbox';
import Grid from '@material-ui/core/Grid';
import Button from '@material-ui/core/Button';
import Paper from '@material-ui/core/Paper';
import Typography from '@material-ui/core/Typography';
import AddIcon from '@material-ui/icons/Add';

import {makeStyles} from '@material-ui/core/styles';

import {QueryGraphContext} from '../../contexts/query/query_graph_context';
import {QueryWhereContext} from '../../contexts/query/query_where_context';
import QueryFilterControl from './query_filter_control';
import {
  QueryFilter,
  QuerySlice,
  EmptyQueryClause
} from '../../contexts/query/query_types';
import QueryClause from './query_clause';
import QueryAnyEverySelectorList from './query_any_every_selector_list';

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
    <QueryClause title='Where'>
      {recordFilters.length > 0 ? (
        <Grid container justify='center' className={classes.header}>
          <Grid item xs={1}>
            OR
          </Grid>
          <Grid item xs={1}>
            Any / Every
          </Grid>
          <Grid item container xs={10}>
            <Grid item xs={2}>
              Model
            </Grid>
            <Grid item xs={9} className={classes.clauseHeading}>
              <Typography>Clauses</Typography>
              <Grid container item>
                <Grid
                  item
                  container
                  spacing={1}
                  alignItems='center'
                  className={classes.clauseSubheading}
                >
                  <Grid item xs={1}>
                    Any/Every
                  </Grid>
                  <Grid item container xs={11}>
                    <Grid item xs={3}>
                      Model
                    </Grid>
                    <Grid item xs={8} container>
                      <Grid item xs={4}>
                        Attribute
                      </Grid>
                      <Grid item xs={3} className={classes.shimLeft}>
                        Operator
                      </Grid>
                      <Grid item xs={4} className={classes.shimLeft}>
                        Operand
                      </Grid>
                      <Grid item xs={1} />
                    </Grid>
                    <Grid item xs={1} />
                  </Grid>
                </Grid>
              </Grid>
            </Grid>
            <Grid item xs={1} />
          </Grid>
        </Grid>
      ) : null}
      {recordFilters.map((filter: QueryFilter, index: number) => (
        <Paper key={index}>
          <Grid
            container
            alignItems='center'
            justify='center'
            className='query-where-selector'
          >
            <Grid item xs={1}>
              <Checkbox
                checked={orRecordFilterIndices.includes(index)}
                color='primary'
                onChange={(e, checked) => handleChangeOrFilters(index)}
                inputProps={{'aria-label': 'secondary checkbox'}}
              />
            </Grid>
            <Grid item xs={1}>
              <QueryAnyEverySelectorList
                filter={filter}
                index={index}
                patchRecordFilter={patchRecordFilter}
              />
            </Grid>
            <Grid item container xs={10} direction='column'>
              <Grid item container alignItems='center' justify='flex-start'>
                <QueryFilterControl
                  key={`${index}-${updateCounter}`}
                  filter={filter}
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
              </Grid>
            </Grid>
          </Grid>
        </Paper>
      ))}
      <Button onClick={addNewRecordFilter} startIcon={<AddIcon />}>
        Filter
      </Button>
    </QueryClause>
  );
};

export default QueryWherePane;
