import React, { useState, useEffect } from 'react';

const ListItem = ({item, onClick }) => {
  let className = 'delete_link';

  // Clone item, but could be simple Object or string
  let display_item = JSON.parse(JSON.stringify(item));

  if (display_item == null || display_item == '') {
    display_item = 'null';
    className = 'delete_link empty';
  } else if (typeof display_item === 'object' &&
    display_item.hasOwnProperty('original_filename')) {
    display_item = display_item.original_filename;
  }

  return(
    <div className='list_item'>
      <span className={ className } onClick={ onClick } >
        { display_item  }
      </span>
    </div>
  );
}

// This is an input to create and edit a list of items
const ListInput = ({ values, itemInput, onChange, onAll, onClear, ...inputProps }) => {
  let [ editing, setEditing ] = useState(null);

  let ItemInput = itemInput;

  const removeValue = (pos) => onChange(
    values.slice(0,pos).concat(values.slice(pos+1))
  );

  const editValue = (new_value) => {
    if (new_value === null || new_value === undefined || new_value === '') return;

    let new_values = values.slice();

    new_values.splice(values.length-1,1,new_value);

    onChange(new_values);
  }

  const addListItem = () => {
    // add a new value to the list
    onChange(values.concat(''));

    // turn on editing
    setEditing(true);
  }

  const addAllItems = () => {
    onAll();

    setEditing(false);
  }

  const removeItems = () => {
    onClear();

    setEditing(false);
  }

  return(
    <div className='list_input'>
      {
        (editing ? values.slice(0,-1) : values).map(
          (value,i) => <ListItem
            key={i}
            item={value}
            onClick={ () => removeValue(i) }
          />
        )
      }
      { editing &&
        <div className='list_item'>
          <ItemInput
            key={values.length}
            onChange={ editValue }
            onBlur={ () => setEditing(null) }
            defaultValue={ values.slice(-1) }
            { ...inputProps }
          />
        </div>
      }
      <div className='list_item'>
        <span className='add_item' onClick={ addListItem }>+</span>
        { onAll ?
          <span className='add_item' onClick={addAllItems}>all</span> : null
        }
        { onClear ?
          <span className='add_item' onClick={removeItems}>reset</span> : null
        }
      </div>
    </div>
  );
}

export default ListInput;
