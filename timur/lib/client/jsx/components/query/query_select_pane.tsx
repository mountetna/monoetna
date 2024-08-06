import React, {useState, useMemo, useContext, useCallback} from 'react';
import Typography from '@material-ui/core/Typography';
import Grid from '@material-ui/core/Grid';
import IconButton from '@material-ui/core/IconButton';
import Tooltip from '@material-ui/core/Tooltip';
import AddIcon from '@material-ui/icons/Add';

import {makeStyles} from '@material-ui/core/styles';

import {DragDropContext, Droppable} from 'react-beautiful-dnd';

import {QueryColumn} from '../../contexts/query/query_types';
import {QueryGraphContext} from '../../contexts/query/query_graph_context';
import {QueryColumnContext} from '../../contexts/query/query_column_context';
import QueryModelAttributeSelector from './query_model_attribute_selector';
import DraggableQueryModelAttributeSelector from './draggable_query_model_attribute_selector';
import QueryClause from './query_clause';
import QueryClauseSummaryControls from './query_clause_summary_controls';
import QueryChevron from './query_chevron';

const useStyles = makeStyles((theme) => ({
  folded: {
    fontStyle: 'italic',
    paddingLeft: '10px',
    cursor: 'pointer'
  },
  columns: {
    width: '100%',
    paddingLeft: '30px'
  }
}));

const QuerySelectPane = () => {
  const {
    state: {graph, rootModel}
  } = useContext(QueryGraphContext);
  const {
    state: {columns},
    addQueryColumn,
    patchQueryColumn,
    removeQueryColumn,
    removeAllQueryColumns,
    setQueryColumns
  } = useContext(QueryColumnContext);
  const classes = useStyles();

  const handleOnSelectModel = useCallback(
    (columnIndex: number, modelName: string, displayLabel: string) => {
      patchQueryColumn(columnIndex, {
        model_name: modelName,
        slices: [],
        attribute_name: '',
        display_label: displayLabel
      });
    },
    [patchQueryColumn]
  );

  const handleOnSelectAttribute = useCallback(
    (columnIndex: number, column: QueryColumn, attributeName: string) => {
      const previousDefaultLabel = `${column.model_name}.${column.attribute_name}`;

      const newLabel = ['', previousDefaultLabel].includes(column.display_label)
        ? `${column.model_name}.${attributeName}`
        : column.display_label;

      patchQueryColumn(columnIndex, {
        model_name: column.model_name,
        slices: [],
        attribute_name: attributeName,
        display_label: newLabel
      });
    },
    [patchQueryColumn]
  );

  const handleOnSelectPredicate = useCallback(
    (columnIndex: number, column: QueryColumn, predicate: string) => {
      patchQueryColumn(columnIndex, {
        ...column,
        predicate
      });
    },
    [patchQueryColumn]
  );

  const handleOnChangeLabel = useCallback(
    (columnIndex: number, column: QueryColumn, label: string) => {
      patchQueryColumn(columnIndex, {
        ...column,
        display_label: label
      });
    },
    [patchQueryColumn]
  );

  const handleOnRemoveColumn = useCallback(
    (columnIndex: number) => {
      removeQueryColumn(columnIndex);
    },
    [removeQueryColumn]
  );

  const modelChoiceSet = useMemo(
    () => [
      ...new Set(
        graph
          .allPaths(rootModel)
          .flat()
          .concat(rootModel ? [rootModel] : [])
      )
    ],
    [graph, rootModel]
  );

  const reorder = (
    columns: QueryColumn[],
    startIndex: number,
    endIndex: number
  ) => {
    const result = Array.from(columns);
    const [removed] = result.splice(startIndex, 1);
    result.splice(endIndex, 0, removed);

    return result;
  };

  const handleOnDragEnd = useCallback(
    (result: any) => {
      if (!result.destination) {
        return;
      }

      if (result.destination.index === result.source.index) {
        return;
      }

      const newColumns = reorder(
        columns,
        result.source.index,
        result.destination.index
      );

      setQueryColumns(newColumns);
    },
    [columns, setQueryColumns]
  );

  const [ fold, setFold ] = useState(true);

  if (!rootModel) return null;

  return (
    <DragDropContext onDragEnd={handleOnDragEnd}>
      <QueryClause title=''>
        <Grid container alignItems='center'>
          <QueryChevron fold={fold} setFold={setFold}/>
          Select attributes as&nbsp;<b>columns</b>:
          <QueryClauseSummaryControls
            fold={fold}
            setFold={setFold}
            addHandler={
              () => addQueryColumn({
                slices: [],
                model_name: '',
                attribute_name: '',
                display_label: ''
              })
            }
            removeHandler={removeAllQueryColumns}
            itemName='column'
            numItems={columns.length}/>
        </Grid>
        { !fold && <Grid container direction='column' className={classes.columns}>
            <Droppable droppableId='columns'>
              {(provided: any) => (
                <div ref={provided.innerRef} {...provided.droppableProps}>
                  {columns.map((column: QueryColumn, index: number) => {
                    const ColumnComponent =
                      0 === index
                        ? QueryModelAttributeSelector
                        : DraggableQueryModelAttributeSelector;
                    return (
                      <ColumnComponent
                        key={index}
                        label='Join Model'
                        column={column}
                        modelChoiceSet={modelChoiceSet}
                        columnIndex={index}
                        canEdit={0 !== index}
                        graph={graph}
                        onSelectModel={(modelName: string) =>
                          handleOnSelectModel(index, modelName, '')
                        }
                        onSelectAttribute={(modelName:string, attributeName: string) =>
                          handleOnSelectAttribute(index, column, modelName, attributeName)
                        }
                        onSelectPredicate={(predicate: string) => {
                          handleOnSelectPredicate(index, column, predicate);
                        }}
                        onChangeLabel={(label: string) =>
                          handleOnChangeLabel(index, column, label)
                        }
                        onRemoveColumn={() => handleOnRemoveColumn(index)}
                        onCopyColumn={() => addQueryColumn({...column})}
                      />
                    );
                  })}
                </div>
              )}
            </Droppable>
          </Grid>
        }
      </QueryClause>
    </DragDropContext>
  );
};

export default QuerySelectPane;
