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
import {Attribute, Model} from 'etna-js/models/magma-model';

const useStyles = makeStyles((theme) => ({
  title: {
    fontWeight: 'bold'
  },
  summary: {
    fontStyle: 'italic',
    paddingLeft: '10px',
  }
}));
const QueryRowPane = () => {
  const {
    state: {graph, rootModel},
    setRootModel
  } = useContext(QueryGraphContext);
  const {setQueryColumns} = useContext(QueryColumnContext);
  const {setDataAndNumRecords} = useContext(QueryResultsContext);

  const onRootModelSelect = useCallback(
    (modelName: string) => {
      let model: Model = graph.models.model(modelName) as Model;
      let column_attrs = sortAttributeList(
        model.selectAttributes(
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
          attribute_name: model.identifier,
          display_label: `${modelName}.${model.identifier}`,
          slices: []
        },
        ...column_attrs.map(
          (attribute:Attribute) => {

          if (isLinkForeignKey(attribute)) {
            const linkModel = graph.models.model(attribute.link_model_name as string) as Model;

            return {
              model_name: attribute.link_model_name,
              attribute_name: linkModel.identifier,
              display_label: `${attribute.link_model_name}.${linkModel.identifier}`,
              slices: []
            }
          }

          return {
            model_name: modelName,
            attribute_name: attribute.attribute_name,
            display_label: attribute.attribute_name,
            slices: []
          }
        })
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
            modelNames={graph.models.modelNames}
            modelName={rootModel || ''}
          />
        </Grid>
      </Grid>
    </QueryClause>
  );
};

export default QueryRowPane;
