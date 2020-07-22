import React, {useMemo, useEffect, useState, useCallback} from 'react';
import {connect} from 'react-redux';

import 'regenerator-runtime/runtime';
import * as _ from 'lodash';

import {css} from '@emotion/core';
import ClimbingBoxLoader from 'react-spinners/ClimbingBoxLoader';

import {
  selectModelNames,
  selectTemplate,
  displayAttributes
} from '../../selectors/magma';
import {
  requestTSV,
  requestModels,
  requestDocuments
} from '../../actions/magma_actions';
import {
  selectSearchCache,
  selectSearchAttributeNames,
  constructSingleFilterString,
  selectSearchFilterParams,
  selectSearchFilterString
} from '../../selectors/search';
import {
  cacheSearchPage,
  setSearchPageSize,
  setSearchPage,
  emptySearchCache,
  setSearchAttributeNames,
  setFilterString
} from '../../actions/search_actions';

import ModelViewer from '../model_viewer';
import useAsyncWork from "etna-js/hooks/useAsyncWork";
import SearchQuery from "./search_query";
import {Loading} from "etna-js/components/Loading";

const spinnerCss = css`
  display: block;
  margin: 2rem auto;
`;

const loadingSpinner =
  <ClimbingBoxLoader
    css={spinnerCss}
    color='green'
    size={20}
    loading={true}
  />


export function Search({
  attribute_names, cache, current_filter, requestDocuments, setSearchPageSize, cacheSearchPage, setSearchPage,
  requestModels, emptySearchCache, setSearchAttributeNames, model_names, requestTSV, magma_state
}) {
  const [pageSize, setPageSize] = useState(10);
  const [selectedModel, setSelectedModel] = useState(undefined);
  const [results, setResults] = useState(0);
  const { current_page, page_size, model_name, record_names, } = cache;
  let { cached_attribute_names } = cache;

  // On mount, essentially.
  useEffect(() => {
    requestModels();
    emptySearchCache();
    setSearchAttributeNames('all');
  }, [])

  const onSelectTableChange = useCallback((model_name) => {
    setSearchAttributeNames('all');
    emptySearchCache();
    setSelectedModel(model_name);
  }, [setSearchAttributeNames, emptySearchCache, setSelectedModel]);

  const [loading, loadDocuments] = useAsyncWork(function* loadDocuments(page, newSearch) {
    const payload = yield requestDocuments({
      model_name: selectedModel,
      record_names: 'all',
      attribute_names: attribute_names,
      filter: current_filter,
      page: page,
      page_size: pageSize,
      collapse_tables: true,
      exchange_name: `request-${selectedModel}`
    });

    if (newSearch) emptySearchCache();

    let model = payload.models[selectedModel];
    setResults(model.count);
    if (!newSearch) setSearchPageSize(pageSize);
    cacheSearchPage(
      page,
      selectedModel,
      Object.keys(model.documents),
      attribute_names,
      page === 1 // clears the cache if you return to page 1
    );
    // Cancel only when multiple consecutive invocations are run.
  }, { cancelWhenChange: [] });

  let pages = model_name ? Math.ceil(results / page_size) : -1;

  const display_attributes = useMemo(() => {
    if (!selectedModel) return [];

    // Have to use the selector here instead of in connect()
    //   because the selectedModel is in component state instead
    //   of global state.
    const template = selectTemplate({ magma: magma_state }, selectedModel);
    return displayAttributes(template);
  }, [selectedModel, magma_state]);

  // We should attempt to re-order the ModelViewer's cached_attribute_names
  //    in the same order as the template's display_attribute_options.
  // This will change in the future once we finalize what a global
  //    attribute ordering should be.
  cached_attribute_names = cached_attribute_names
    ? _.flatten(
      display_attributes.filter(
        (opt) =>
          // This will certainly change when it's an array of items
          cached_attribute_names.includes(opt) ||
          cached_attribute_names === 'all'
      )
    )
    : null;

  return (
    <div id='search'>
      <div className='control'>
        <SearchQuery loading={loading} onSelectTableChange={onSelectTableChange} pageSize={pageSize}
                     display_attributes={display_attributes}
                     selectedModel={selectedModel} setPage={setPage} setPageSize={setPageSize} />
        <Loading loading={results === 0 || loading}>
          <div className='results'>
            Found {results} records in{' '}
            <span className='model_name'>{model_name}</span>
          </div>
        </Loading>
      </div>
      <div className='body'>
        <Loading loading={!model_name || (loading && loadingSpinner)} delay={400}>
          <div className='documents'>
            <ModelViewer
              model_name={model_name}
              record_names={record_names}
              page={current_page - 1}
              pages={pages}
              page_size={page_size}
              setPage={setPage}
              restricted_attribute_names={
                cached_attribute_names !== 'all'
                  ? cached_attribute_names
                  : null
              }
            />
          </div>
        </Loading>
      </div>
    </div>
  );

  function setPage(page, newSearch) {
    // The page model offset + 1
    page++;
    // Need to re-fetch a page if the user has clicked a new set of
    //    attribute names from the TreeView
    if (!cache.isCached(page.toString()) || newSearch) {
      loadDocuments(page, newSearch)
    }
    setSearchPage(page)
  }
}

export default connect(
  (state, props) => ({
    model_names: selectModelNames(state),
    cache: selectSearchCache(state),
    attribute_names: selectSearchAttributeNames(state),
    current_filter: constructSingleFilterString(state),
    filter_string: selectSearchFilterString(state),
    filter_params: selectSearchFilterParams(state),
    magma_state: state.magma
  }),
  {
    requestModels,
    cacheSearchPage,
    setSearchPage,
    setSearchPageSize,
    setSearchAttributeNames,
    setFilterString,
    emptySearchCache,
    requestDocuments,
    requestTSV
  }
)(Search);
