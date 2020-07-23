import React, {useState} from 'react';
import {connect} from 'react-redux';
import SelectInput from '../inputs/select_input';
import {selectModelNames} from "../../selectors/magma";
import {
  constructSingleFilterString,
  selectSearchAttributeNames,
} from "../../selectors/search";
import {requestTSV} from "../../actions/magma_actions";
import {
  setFilterString, setSearchAttributeNames,
} from "../../actions/search_actions";
import CollapsibleArrow from "etna-js/components/CollapsibleArrow";
import QueryBuilder from "./query_builder";

export function SearchQuery({
  selectedModel, setFilterString, loading, requestTSV, model_names, onSelectTableChange,
  pageSize, setPageSize, setPage, current_filter, attribute_names, display_attributes,
}) {
  const buttonDisabled = !selectedModel || loading;
  const buttonClasses = buttonDisabled ? 'button disabled' : 'button';
  const [showAdvanced, setShowAdvanced] = useState(false);

  const tableOptionsLine = <React.Fragment>
    <span className='label'>Show table</span>
    <SelectInput
      name='model'
      values={model_names}
      value={selectedModel}
      onChange={onSelectTableChange}
      showNone='enabled'
    />

    <span className='label'>Page size</span>
    <SelectInput
      values={[10, 25, 50, 200]}
      defaultValue={pageSize}
      onChange={setPageSize}
      showNone='disabled'
    />

    <input
      id='search-pg-search-btn'
      type='button'
      className={buttonClasses}
      value='Search'
      disabled={buttonDisabled}
      onClick={() => setPage(0, true)}
    />
    <input
      id='search-pg-tsv-btn'
      className={buttonClasses}
      type='button'
      value={'\u21af TSV'}
      disabled={buttonDisabled}
      onClick={() =>
        requestTSV(
          selectedModel,
          current_filter,
          attribute_names
        )
      }
    />
  </React.Fragment>;

  const advancedSearch = <input
    type='text'
    className='filter'
    placeholder='Filter query'
    onBlur={(e) => setFilterString(e.target.value)}
  />;

  return (
    <div className='query'>
      <CollapsibleArrow label={tableOptionsLine}>
        { selectedModel && <div>
          {showAdvanced ? advancedSearch : <QueryBuilder selectedModel={selectedModel} display_attributes={display_attributes} /> }
        </div> }
      </CollapsibleArrow>
    </div>
  );
}


export default connect(
  (state) => ({
    model_names: selectModelNames(state),
    attribute_names: selectSearchAttributeNames(state),
    current_filter: constructSingleFilterString(state),
  }),
  {
    setSearchAttributeNames,
    setFilterString,
    requestTSV,
  }
)(SearchQuery);

