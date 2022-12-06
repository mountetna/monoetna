import React, {useContext} from 'react';

import {useReduxState} from 'etna-js/hooks/useReduxState';
import QueryControls from './query_controls';
import QueryResults from './query_results';
import {QueryGraphContext} from '../../contexts/query/query_graph_context';

import useQueryGraph from '../../contexts/query/use_query_graph';

const QueryBuilder = ({
  syncQueryParams,
  queryControlButtons
}: {
  syncQueryParams: boolean;
  queryControlButtons: React.FC;
}) => {
  const {
    state: {graph},
    setGraph
  } = useContext(QueryGraphContext);
  const reduxState = useReduxState();

  useQueryGraph(reduxState, graph, setGraph);

  if (!graph || !graph.initialized) return null;

  return (
    <React.Fragment>
      <QueryControls
        syncQueryParams={syncQueryParams}
        queryControlButtons={queryControlButtons}
      />
      <QueryResults />
    </React.Fragment>
  );
};

export default QueryBuilder;
