import React, {useState, useCallback, useEffect} from 'react';

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

  return (
    <div className='add-attribute-modal'>
      <div className='header'>Add Attribute</div>
      <div className='options-tray tray'>
        <label>
          <input onChange={setName} type='text' />
          Attribute name
        </label>
        <label>
          <input onChange={setDescription} type='text' />
          Attribute description
        </label>
        <label>
          <input onChange={setType} type='text' />
          Attribute type
        </label>
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
