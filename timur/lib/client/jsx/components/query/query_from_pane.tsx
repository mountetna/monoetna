import React, {useCallback, useContext} from 'react';

import Grid from '@material-ui/core/Grid';
import IconButton from '@material-ui/core/IconButton';

import {QueryGraphContext} from '../../contexts/query/query_graph_context';
import {QueryColumnContext} from '../../contexts/query/query_column_context';
import {QueryResultsContext} from '../../contexts/query/query_results_context';
import {EmptyQueryResponse} from '../../contexts/query/query_types';
import QueryModelSelector from './query_model_selector';
import QueryClause from './query_clause';
import QueryChevron from './query_chevron';
import QueryControlButtons from './query_control_buttons';
import {isLinkCollection, sortAttributeList} from '../../utils/attributes';

const QueryFromPane = () => {
  const {
    state: {graph, rootModel},
    setRootModel
  } = useContext(QueryGraphContext);
  const {setQueryColumns} = useContext(QueryColumnContext);
  const {setDataAndNumRecords} = useContext(QueryResultsContext);

  const onRootModelSelect = useCallback(
    (modelName: string) => {
      let template = graph.template(modelName);
      let column_attrs = sortAttributeList(
        Object.values(template.attributes).filter(
          attribute => !(isLinkCollection( attribute )
            || attribute.hidden
            || attribute.attribute_type == 'identifier'
          )
        ),
        true
      );
      setRootModel(modelName);
      setQueryColumns([
        {
          model_name: modelName,
          attribute_name: template.identifier,
          display_label: `${modelName}.${template.identifier}`,
          slices: []
        },
        ...column_attrs.map(
          attribute => ({
            model_name: modelName,
            attribute_name: attribute.attribute_name,
            display_label: attribute.attribute_name,
            slices: []
          })
        )
      ]);
      setDataAndNumRecords(EmptyQueryResponse, 0);
    },
    [graph, setRootModel, setQueryColumns, setDataAndNumRecords]
  );

  return (
    <QueryClause title='From'>
      <Grid item container>
        <Grid item container xs={8} alignItems='center'>
          <QueryChevron disabled/>
          Select records from model&nbsp;
          <QueryModelSelector
            setModel={onRootModelSelect}
            modelNames={[...graph.allowedModels]}
            modelName={rootModel || ''}
          />&nbsp;as&nbsp;<b>rows</b>
        </Grid>
        <Grid item container alignItems='center' justify='flex-end' xs={4}>
          {rootModel ? <QueryControlButtons /> : null}
        </Grid>
      </Grid>
    </QueryClause>
  );
};

export default QueryFromPane;
