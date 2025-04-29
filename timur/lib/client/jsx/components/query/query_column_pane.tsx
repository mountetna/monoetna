import React, {useEffect, useState, useMemo, useContext, useCallback} from 'react';
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
import QueryColumnSelector from './query_column_selector';
import DraggableQueryColumnSelector from './draggable_query_column_selector';
import QueryClause from './query_clause';
import QueryClauseSummaryControls from './query_clause_summary_controls';
import QueryChevron from './query_chevron';
import MapSelector from './map_selector';
import {isLink} from '../../utils/attributes';
import {Attribute} from 'etna-js/models/magma-model';

const useStyles = makeStyles((theme) => ({
  folded: {
    fontStyle: 'italic',
    paddingLeft: '10px',
    cursor: 'pointer'
  },
  title: {
    fontWeight: 'bold'
  },
  columns: {
    width: '100%',
    paddingLeft: '30px'
  }
}));

const QueryColumnPane = () => {
  const { state: {graph, rootModel} } = useContext(QueryGraphContext);
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
    (columnIndex: number, column: QueryColumn, modelName: string, attributeName: string) => {
      const previousDefaultLabel = `${column.model_name}.${column.attribute_name}`;

      const newLabel = ['', previousDefaultLabel].includes(column.display_label)
        ? `${modelName}.${attributeName}`
        : column.display_label;

      patchQueryColumn(columnIndex, {
        model_name: modelName,
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
      ...graph.connectedModelsAnd(rootModel)
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

      setQueryColumns(newColumns as QueryColumn[]);
    },
    [columns, setQueryColumns]
  );

  const [ columnsModel, setColumnsModel ] = useState(rootModel);

  useEffect(() => {
    setColumnsModel(rootModel);
  }, [ rootModel ]);

  const handleOnSelectAttributes = useCallback(
    (attribute_names: string[]) => {
      const newColumns = [ ...columns, ...attribute_names.map(
        attribute_name => ({
          slices: [],
          model_name: columnsModel,
          attribute_name,
          display_label: attribute_name
        })
      ) ];
      setQueryColumns(newColumns as QueryColumn[]);
    },
    [columns, setQueryColumns, columnsModel]
  );
    
  const [ fold, setFold ] = useState(true);

  const [ showAttributesModal, setShowAttributesModal ] = useState(false);

  if (!rootModel) return null;

  return (
    <DragDropContext onDragEnd={handleOnDragEnd}>
      <QueryClause title=''>
        <Grid container alignItems='center'>
          <QueryChevron fold={fold} setFold={setFold}/>
          <Typography className={classes.title}>Columns:</Typography>
          <QueryClauseSummaryControls
            fold={fold}
            setFold={setFold}
            addHandler={ () => setShowAttributesModal(true) }
            removeHandler={removeAllQueryColumns}
            itemName='column'
            numItems={columns.length}/>
          <MapSelector
            open={showAttributesModal}
            onClose={() => setShowAttributesModal(false)}
            setAttributes={handleOnSelectAttributes}
            modelName={columnsModel as string}
            modelNames={modelChoiceSet}
            filterAttributes={
              (attribute:Attribute) => !(isLink(attribute) || attribute.hidden)
            }
            setModel={setColumnsModel}
          />
        </Grid>
        { !fold && <Grid container direction='column' className={classes.columns}>
            <Droppable droppableId='columns'>
              {(provided: any) => (
                <div ref={provided.innerRef} {...provided.droppableProps}>
                  {columns.map((column: QueryColumn, index: number) => {
                    const ColumnComponent =
                      0 === index
                        ? QueryColumnSelector
                        : DraggableQueryColumnSelector;
                    return (
                      <ColumnComponent
                        key={index}
                        label='Join Model'
                        column={column}
                        modelChoiceSet={modelChoiceSet}
                        columnIndex={index}
                        canEdit={0 !== index}
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

export default QueryColumnPane;
