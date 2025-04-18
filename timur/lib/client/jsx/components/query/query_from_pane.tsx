import React, {useCallback, useContext} from 'react';

import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
import IconButton from '@material-ui/core/IconButton';
import {makeStyles} from '@material-ui/core/styles';

import {QueryGraphContext} from '../../contexts/query/query_graph_context';
import {QueryColumnContext} from '../../contexts/query/query_column_context';
import {QueryResultsContext} from '../../contexts/query/query_results_context';
import {EmptyQueryResponse} from '../../contexts/query/query_types';
import QueryModelSelector from './query_model_selector';
import QueryClause from './query_clause';
import QueryChevron from './query_chevron';
import {isLinkForeignKey, isLinkCollection, sortAttributeList} from '../../utils/attributes';
import {Attribute, Template} from 'etna-js/models/magma-model';

const useStyles = makeStyles((theme) => ({
  title: {
    fontWeight: 'bold'
  },
  summary: {
    fontStyle: 'italic',
    paddingLeft: '10px',
  }
}));
const QueryFromPane = () => {
  const {
    state: {graph, rootModel},
    setRootModel
  } = useContext(QueryGraphContext);
  const {setQueryColumns} = useContext(QueryColumnContext);
  const {setDataAndNumRecords} = useContext(QueryResultsContext);

  const onRootModelSelect = useCallback(
    (modelName: string) => {
      let template: Template = graph.template(modelName);
      let column_attrs = sortAttributeList(
        Object.values(template.attributes).filter(
          (attribute: Attribute) => !(isLinkCollection( attribute )
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
          (attribute:Attribute) => (isLinkForeignKey(attribute) ? {
            model_name: attribute.link_model_name,
            attribute_name: graph.template(attribute.link_model_name as string).identifier,
            display_label: `${attribute.link_model_name}.${graph.template(attribute.link_model_name as string).identifier}`,
            slices: []
          } : {
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

  const classes = useStyles();

  return (
    <QueryClause title=''>
      <Grid item container>
        <Grid item container xs={8} alignItems='center'>
          <QueryChevron disabled/>
          <Typography className={classes.title}>Rows:</Typography>
          <Typography className={classes.summary}>from&nbsp;</Typography>
          <QueryModelSelector
            setModel={onRootModelSelect}
            modelNames={Object.keys(graph.models)}
            modelName={rootModel || ''}
          />
        </Grid>
      </Grid>
    </QueryClause>
  );
};

export default QueryFromPane;
