import React, {useState, useCallback, useEffect} from 'react';

import TextField from '@material-ui/core/TextField';
import MenuItem from '@material-ui/core/MenuItem';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {useModal} from 'etna-js/components/ModalDialogContainer';

import DisabledButton from '../search/disabled_button';

export default function AddAttributeModal({onSave}: {onSave: any}) {
  const [disabled, setDisabled] = useState(true);
  const [name, setName] = useState(null);
  const [description, setDescription] = useState(null);
  const [type, setType] = useState(null);
  const {dismissModal} = useModal();
  const invoke = useActionInvoker();

  const handleOnSave = useCallback(() => {
    onSave({
      name,
      description,
      type
    });
  }, [name, description, type]);

  useEffect(() => {
    if (name && description && type) {
      setDisabled(false);
    } else {
      setDisabled(true);
    }
  }, [name, description, type]);

  const handleOnCancel = useCallback(() => {
    invoke(dismissModal());
  }, [invoke, dismissModal]);

  const attributeTypes = ['string', 'shifted_date_time', 'date_time'];

  return (
    <div className='add-attribute-modal'>
      <div className='header'>Add Attribute</div>
      <div className='options-tray tray'>
        <TextField id='attribute-name' label='Name' onChange={setName} />
        <TextField
          id='attribute-description'
          label='Description'
          onChange={setDescription}
        />
        <TextField
          id='attribute-type'
          select
          value=''
          label='Type'
          onChange={(e: any) => setType(e.target.value)}
        >
          {attributeTypes.sort().map((option, i) => (
            <MenuItem key={i} value={option}>
              {option}
            </MenuItem>
          ))}
        </TextField>
      </div>
      <div className='options-action-wrapper'>
        <DisabledButton
          id='cancel-add-attribute-btn'
          className='cancel'
          label='Cancel'
          disabled={false}
          onClick={handleOnCancel}
        />
        <DisabledButton
          id='add-attribute-btn'
          className='save'
          label='Save'
          disabled={disabled}
          onClick={handleOnSave}
        />
      </div>
    </div>
  );
}
