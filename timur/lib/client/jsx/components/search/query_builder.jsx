import React, {useCallback, useState, useEffect} from 'react';
import {connect} from "react-redux";
import {selectSearchAttributeNames} from "../../selectors/search";
import {setFilterString, setSearchAttributeNames} from "../../actions/search_actions";
import {useModal} from "etna-js/components/ModalDialogContainer";
import TreeView, {getSelectedLeaves} from 'etna-js/components/TreeView';
import SelectInput from "../inputs/select_input";

export function QueryBuilder({ display_attributes, setFilterString }) {
  const { openModal } = useModal();
  const [filtersState, setFiltersState] = useState([]);

  useEffect(() => {
    // TODO: Compse the filter string from filtersState and calling setFilterString
  }, [setFilterString, filtersState]);

  const onOpenAttributeFilter = () => {
    openModal(<FilterAttributesModal display_attributes={display_attributes}/>);
  };

  const onOpenFilters = () => {
    openModal(<QueryFilterModal filtersState={filtersState} setFiltersState={setFilterString} />);
  }

  return <div className='query-builder'>
    <a className='pointer' onClick={onOpenAttributeFilter}>
      Add/Remove Attributes
    </a>
    <a className='pointer' onClick={onOpenFilters}>
      Add/Remove Filters
    </a>
  </div>;
}

function FilterAttributesModal({ setSearchAttributeNames, display_attributes, attribute_names }) {
  const display_attribute_options = [['All', display_attributes.map(e => [e])]];
  const { dismissModal } = useModal();
  const handleTreeViewSelectionsChange = useCallback((new_state) => {
    setSearchAttributeNames(getSelectedLeaves(new_state));
  }, [setSearchAttributeNames]);


  return (
    <div className='search-attribute-filters-modal'>
      <TreeView
        flowHorizontally={true}
        selected={attributeNamesToSelected(attribute_names)}
        options={display_attribute_options}
        onSelectionsChange={handleTreeViewSelectionsChange}
      />
      <div className='actions'>
        <button onClick={dismissModal} disabled={attribute_names.length === 0}>Ok</button>
      </div>
    </div>
  );
}

FilterAttributesModal = connect(
  (state) => ({ attribute_names: selectSearchAttributeNames(state), }),
  { setSearchAttributeNames, }
)(FilterAttributesModal);

const operators = [
  '>',
  '<'
];

function QueryFilterModal({ attribute_names, filterState: initialFilterState, setFilterState: updateParentFilterState }) {
  const { dismissModal } = useModal();
  let [filterState, setLocalFilterState] = useState(initialFilterState);

  function setFilterState(state) {
    setLocalFilterState(state);
    updateParentFilterState(state);
  }

  if (attribute_names === 'all') {
    // TODO HERE
    attribute_names = [];
  }

  // TODO HERE
    function onFilterAttributeChange(idx) {
      return function (value) {
        if (!value) {
          filterState = [...filterState];
          filterState.splice(idx, 1);
          setFilterState(filterState);
        } else {
          if (idx === -1) {
            filterState = [...filterState, { attribute: value }];
            setFilterState(filterState);
          } else {
            filterState = [...filterState];
            filterState.splice(idx, 1, { ...filterState[idx], attribute: value });
          }
        }
      }
    }

  return (
    <div className='search-filters-modal'>
      <div>
        <SelectInput
          name='model'
          values={attribute_names}
          showNone='enabled'
        />
        <SelectInput
          name='model'
          values={operators}
          showNone='enabled'
        />
      </div>
      <div className='actions'>
        <button onClick={dismissModal} disabled={attribute_names.length === 0}>Ok</button>
      </div>
    </div>
  );
}

QueryFilterModal = connect(
  (state) => ({ attribute_names: selectSearchAttributeNames(state), }),
)(QueryFilterModal);

function attributeNamesToSelected(attributeNames) {
  if (attributeNames === 'all') return null;
  return { All: attributeNames.reduce((o, p) => (o[p] = true, o), {}) };
}

export default connect(
  (state) => ({}),
  {
    setFilterString,
  }
)(QueryBuilder);
