import React, {useMemo} from 'react';
import _ from 'lodash';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {selectAllowedModelAttributes} from '../../selectors/query_selector';
import {visibleSortedAttributesWithUpdatedAt} from '../../utils/attributes';
import {QueryGraph} from '../../utils/query_graph';

const useQueryClause = ({
  modelName,
  graph,
  isColumnFilter
}: {
  modelName: string;
  graph: QueryGraph;
  isColumnFilter: boolean;
}) => {
  const invoke = useActionInvoker();

  const modelAttributes = useMemo(() => {
    if ('' !== modelName) {
      const model = graph.models.model(modelName);
      if (!model) return [];

      let sortedTemplateAttributes = visibleSortedAttributesWithUpdatedAt(
        model.attributes
      );

      return selectAllowedModelAttributes(
        sortedTemplateAttributes,
        !isColumnFilter
      );
    }
    return [];
  }, [modelName, graph, isColumnFilter]);

  return {
    modelAttributes
  };
};

export default useQueryClause;
