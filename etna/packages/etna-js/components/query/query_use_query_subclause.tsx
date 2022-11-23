import React, {useMemo, useState, useCallback} from 'react';
import _ from 'lodash';

import {Cancellable} from 'etna-js/utils/cancellable';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {showMessages} from 'etna-js/actions/message_actions';
import {requestAnswer} from 'etna-js/actions/magma_actions';
import {QuerySubclause} from '../../contexts/query/query_types';
import {QueryGraph} from '../../utils/query_graph';

const useQuerySubclause = ({
  subclause,
  modelName,
  graph
}: {
  subclause: QuerySubclause;
  modelName: string;
  graph: QueryGraph;
}) => {
  const [distinctAttributeValues, setDistinctAttributeValues] = useState(
    [] as string[]
  );
  const invoke = useActionInvoker();

  const attributeType = useMemo(() => {
    if ('' !== subclause.attributeName) {
      const template = graph.template(modelName);
      if (!template) return '';

      return template.attributes[
        subclause.attributeName
      ].attribute_type.toLowerCase();
    }
    return '';
  }, [subclause.attributeName, modelName, graph]);

  const fetchDistinctAttributeValues = useCallback(() => {
    const cancellable = new Cancellable();

    if ('' === modelName || '' === subclause.attributeName) {
      setDistinctAttributeValues([]);
    } else if ('string' !== subclause.attributeType) {
      setDistinctAttributeValues([]);
    } else {
      cancellable
        .race(
          invoke(
            requestAnswer({
              query: [modelName, '::distinct', subclause.attributeName]
            })
          )
        )
        .then(({result, cancelled}: any) => {
          if (result && !cancelled) setDistinctAttributeValues(result.answer);
        })
        .catch((e: any) => {
          invoke(showMessages([e]));
        });
    }

    return () => cancellable.cancel();
  }, [subclause, invoke, modelName]);

  return {
    attributeType,
    fetchDistinctAttributeValues,
    distinctAttributeValues
  };
};

export default useQuerySubclause;
