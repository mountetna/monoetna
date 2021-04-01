import { connect } from 'react-redux';

import React, { useState, useEffect } from 'react';
import { requestPredicates, requestModels } from 'etna-js/actions/magma_actions';
import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import PredicateChainSet from './query_builder/predicate_chain_set';

// helper to format the predicate state into an actual query
const predicateArray = (predicate) => {
  let { type, filters, model_name, args, start } = predicate;
  switch(type) {
    case 'model':
      let ary = [];
      if (start) ary.push(model_name);
      if (filters.length) ary = ary.concat( filters.map(chainArray) );
      return ary.concat(args);
    case 'terminal':
      return [];
    case 'record':
      if (Array.isArray(args[0])) return [ args[0].map(chainArray) ];
      // let this continue through to default
    default:
      return args.filter(x=>x!=null);
  }
};

// helper to format a series of predicates
const chainArray = (chain) => chain.map(predicateArray).reduce(
  (predArray,pred) => predArray.concat(pred),
  []
);

// test whether the predicate is 'done' and we should calculate a query string to display
const predicateComplete = (predicate) => predicate.completed || predicate.type == 'terminal';

// generate an actual string from a predicate array
const formatChainArray = (terms) => {
  return `[ ${
    terms.map(term => {
      if (typeof(term) == 'string') return `'${term}'`;
      if (Array.isArray(term)) return formatChainArray(term);
      return term;
    }).join(', ')
    } ]`;
};

// initially there is an empty model_list predicate
const defaultQuery = () => (
  [
    [
      {
        type: 'model',
        start: true
      }
    ]
  ]
);

// The root component that renders the query builder - this holds the state of the current query
const QueryBuilder = () => {
  let [ query, updateQuery ] = useState(defaultQuery());
  const invoke = useActionInvoker();

  useEffect( () => {
    invoke(requestModels());
    invoke(requestPredicates());
  });

  let predicates = query[0];

  return <div id='query'>
    <PredicateChainSet chains={ query } update={ updateQuery } />
    <div className='query'>
      {
        predicates.every(predicateComplete)
          ? formatChainArray(chainArray(predicates))
          : null
      }
    </div>
  </div>;
}

export default QueryBuilder;
