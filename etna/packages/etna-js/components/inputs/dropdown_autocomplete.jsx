// Wraps the Dropdown component and provides a text
//   input, so user can type into the input field and
//   see a list of filtered options.
import React, {useEffect, useState} from 'react';

import Icon from '../icon';
import debounce from '../../utils/debounce';

export default function DropdownAutocomplete({
  list,
  onSelect,
  defaultValue,
  waitTime
}) {
  var collator = new Intl.Collator(undefined, {
    numeric: true,
    sensitivity: 'base'
  });

  const sortedList = list.sort(collator.compare);

  const [filteredList, setFilteredList] = useState(null);
  const [showList, setShowList] = useState(false);
  const [selectedValue, setSelectedValue] = useState(null);

  function filterTheList(value) {
    let re = new RegExp(value);
    setFilteredList(
      sortedList.filter((item) => item.match(re)).slice(0, maxItems || 10)
    );
  }

  function onSelectItem(value) {
    console.log('on select', value);
    onSelect(value);
    setSelectedValue(value);
    setShowList(false);
  }

  function openList() {
    setShowList(true);
  }

  function closeList() {
    // Clicking on an <li> triggers onBlur
    //   for the <SlowTextInput> first,
    //   so we set a timeout to wait for
    //   that event to happen.
    setTimeout(() => {
      setShowList(false);
    }, 150);
  }

  function toggleList() {
    showList ? closeList() : openList();
  }

  function onChange(value) {
    filterTheList(value);
    setSelectedValue(value);
    setShowList(true);
  }

  function handleChange(e) {
    onChange(e.target.value);
  }

  useEffect(() => {
    onChange = debounce(onChange, waitTime || 500);
    setSelectedValue(defaultValue ? defaultValue : '');
  }, []);

  useEffect(() => {
    setFilteredList(sortedList || []);
  }, [showList]);

  return (
    <div className='dropdown-autocomplete-input'>
      <div className='dropdown-autocomplete-input-field'>
        <input type='text' onChange={handleChange} value={selectedValue} />
        <div className='icon-wrapper' onClick={toggleList}>
          <Icon icon={`${showList ? 'caret-up' : 'caret-down'}`}></Icon>
        </div>
      </div>
      {showList ? (
        <ul className={`dropdown-autocomplete-options`}>
          {filteredList && filteredList.length > 0 ? (
            filteredList.slice(0, maxItems || 10).map((item, index) => (
              <li onClick={() => onSelectItem(item)} key={index}>
                {item}
              </li>
            ))
          ) : (
            <li> --- No results --- </li>
          )}
        </ul>
      ) : null}
    </div>
  );
}
