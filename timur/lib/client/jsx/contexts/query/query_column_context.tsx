import React, {useState, createContext, useCallback} from 'react';

import {QueryColumn} from './query_types';
import {QueryGraph} from '../../utils/query_graph';
import {Attribute, Model, Models} from 'etna-js/models/magma-model';
import {isLinkForeignKey, isLinkCollection, sortAttributeList} from '../../utils/attributes';

export const defaultQueryColumnParams = {
  columns: [] as QueryColumn[]
};

const defaultQueryColumnState = {
  ...defaultQueryColumnParams
};

export type QueryColumnState = Readonly<typeof defaultQueryColumnState>;

export const defaultQueryColumnContext = {
  state: defaultQueryColumnState as QueryColumnState,
  addQueryColumn: (column: QueryColumn) => {},
  removeQueryColumn: (index: number) => {},
  removeAllQueryColumns: (rootModel: string, graph: QueryGraph) => {},
  patchQueryColumn: (index: number, column: QueryColumn) => {},
  setRootIdentifierColumn: (column: QueryColumn) => {},
  setQueryColumns: (columns: QueryColumn[]) => {}
};

export type QueryColumnContextData = typeof defaultQueryColumnContext;

export const QueryColumnContext = createContext(defaultQueryColumnContext);
export type QueryColumnContext = typeof QueryColumnContext;
export type ProviderProps = {params?: {}; children: any};

export const columnsForModel = (model: Model, models: Models) => {
  let column_attrs = sortAttributeList(
    model.selectAttributes(
      (attribute: Attribute) => !(isLinkCollection( attribute )
        || attribute.hidden
        || attribute.attribute_type == 'identifier'
      )
    ),
    true
  )
  return [
    {
      model_name: model.name,
      attribute_name: model.identifier,
      display_label: `${model.name}.${model.identifier}`,
      slices: []
    },
    ...column_attrs.map(
      (attribute:Attribute) => {
        if (isLinkForeignKey(attribute)) {
          const linkModel = models.model(attribute.link_model_name as string) as Model;

          return {
            model_name: attribute.link_model_name,
            attribute_name: linkModel.identifier,
            display_label: `${attribute.link_model_name}.${linkModel.identifier}`,
            slices: []
          }
        }

        return {
          model_name: model.name,
          attribute_name: attribute.attribute_name,
          display_label: attribute.attribute_name,
          slices: []
        }
      }
    )
  ];
}

export const QueryColumnProvider = (
  props: ProviderProps & Partial<QueryColumnContextData>
) => {
  const [state, setState] = useState(
    props.state || defaultQueryColumnContext.state
  );

  const addQueryColumn = useCallback(
    (column: QueryColumn) => {
      setState({
        ...state,
        columns: [...(state.columns || [])].concat([column])
      });
    },
    [state]
  );

  const removeQueryColumn = useCallback(
    (index: number) => {
      let updatedQueryColumns = [...state.columns];
      updatedQueryColumns.splice(index, 1);
      setState({
        ...state,
        columns: updatedQueryColumns
      });
    },
    [state]
  );

  const removeAllQueryColumns = useCallback(
    (rootModel, graph) => {
      setState({
        ...state,
        columns: [
          state.columns[0],
          ...(
            graph.models.model(rootModel).isTable
              ? [ state.columns[1] ]
              : []
          )
        ]
      });
    },
    [state]
  );

  const patchQueryColumn = useCallback(
    (index: number, column: QueryColumn) => {
      let updatedQueryColumns = [...state.columns];
      updatedQueryColumns[index] = column;
      setState({
        ...state,
        columns: updatedQueryColumns
      });
    },
    [state]
  );

  const setRootIdentifierColumn = useCallback(
    (column: QueryColumn) => {
      setState({
        ...state,
        columns: [column]
      });
    },
    [state]
  );

  const setQueryColumns = useCallback(
    (columns: QueryColumn[]) => {
      setState({
        ...state,
        columns
      });
    },
    [state]
  );

  return (
    <QueryColumnContext.Provider
      value={{
        state,
        addQueryColumn,
        removeQueryColumn,
        removeAllQueryColumns,
        patchQueryColumn,
        setRootIdentifierColumn,
        setQueryColumns
      }}
    >
      {props.children}
    </QueryColumnContext.Provider>
  );
};
