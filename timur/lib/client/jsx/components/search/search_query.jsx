import React, {useState} from 'react';
import {connect} from 'react-redux';
import SelectInput from 'etna-js/components/inputs/select_input';
import Toggle from 'etna-js/components/inputs/toggle';
import {selectModelNames} from 'etna-js/selectors/magma';
import {useModal} from 'etna-js/components/ModalDialogContainer';
import {selectSearchFilterString} from '../../selectors/search';
import {
  setFilterString,
  setSearchAttributeNames
} from '../../actions/search_actions';
import QueryBuilder from './query_builder';
import DisabledButton from './disabled_button';
import SearchTsvOptionsModal from './search_tsv_options_modal';

const TableSelect = ({model_names, selectedModel, onSelectTableChange}) => (
  <div className='table-select'>
    <span className='label'>Search table</span>
    <SelectInput
      name='model'
      values={model_names}
      value={selectedModel}
      onChange={onSelectTableChange}
      showNone='enabled'
    />
  </div>
);

const PageSelect = ({pageSize, setPageSize}) => (
  <div className='page-select'>
    <span className='label'>Page size</span>
    <SelectInput
      values={[10, 25, 50, 200]}
      defaultValue={pageSize}
      onChange={setPageSize}
      showNone='disabled'
    />
  </div>
);

export function SearchQuery({
  selectedModel,
  setFilterString,
  loading,
  model_names,
  onSelectTableChange,
  pageSize,
  setPageSize,
  setPage,
  display_attributes,
  filter_string
}) {
  const buttonDisabled = !selectedModel || loading;
  const [showAdvanced, setShowAdvanced] = useState(false);

  const {openModal} = useModal();

  const advancedSearch = (
    <div className='advanced-search'>
      <input
        type='text'
        className='filter'
        placeholder='Filter query'
        defaultValue={filter_string}
        onBlur={(e) => setFilterString(e.target.value)}
      />
    </div>
  );

  return (
    <div className='query'>
      <TableSelect {...{selectedModel, model_names, onSelectTableChange}} />
      <PageSelect {...{pageSize, setPageSize}} />

      <div className='query-options'>
        {selectedModel && (
          <div className='query-mode'>
            <Toggle
              label={showAdvanced ? 'Raw' : 'Basic'}
              selected={showAdvanced}
              onClick={() => setShowAdvanced(!showAdvanced)}
            />
          </div>
        )}
        {selectedModel &&
          (showAdvanced ? (
            advancedSearch
          ) : (
            <QueryBuilder
              setShowAdvanced={setShowAdvanced}
              selectedModel={selectedModel}
              display_attributes={display_attributes}
            />
          ))}
      </div>
      <DisabledButton
        id='search-pg-search-btn'
        label='Search'
        disabled={buttonDisabled}
        onClick={() => setPage(0, true)}
      />
      <DisabledButton
        id='search-pg-tsv-btn'
        label={'\u21af TSV'}
        disabled={buttonDisabled}
        onClick={() =>
          openModal(<SearchTsvOptionsModal selectedModel={selectedModel} />)
        }
      />
    </div>
  );
}

export default connect(
  (state) => ({
    model_names: selectModelNames(state),
    filter_string: selectSearchFilterString(state)
  }),
  {
    setSearchAttributeNames,
    setFilterString
  }
)(SearchQuery);
