import React, {useState, useCallback, useEffect} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {useModal} from 'etna-js/components/ModalDialogContainer';

import DisabledButton from 'etna-js/components/disabled_button';
import {ShrinkingLabelTextField} from './shrinking_label_text_field';
import {COMMA_SEP, SNAKE_CASE} from '../../utils/edit_map';
import ModalSelect from './modal_select';

export default function AddAttributeModal({onSave}: {onSave: any}) {
  const [disabled, setDisabled] = useState(true);
  const [name, setName] = useState('');
  const [description, setDescription] = useState('');
  const [type, setType] = useState('');
  const [group, setGroup] = useState('');
  const {dismissModal} = useModal();
  const invoke = useActionInvoker();

  const handleOnSave = useCallback(() => {
    onSave({
      attribute_name: name,
      description,
      type,
      attribute_group: group
    });
  }, [name, description, type, group]);

  useEffect(() => {
    if (name && type) {
      setDisabled(false);
    } else {
      setDisabled(true);
    }
  }, [name, description, type]);

  const handleOnCancel = useCallback(() => {
    invoke(dismissModal());
  }, [invoke, dismissModal]);

  const attributeTypes = [
    'string',
    'shifted_date_time',
    'date_time',
    'boolean',
    'file',
    'image',
    'file_collection',
    'float',
    'integer'
  ];

  return (
    <div className='add-attribute-modal model-actions-modal'>
      <div className='header'>Add Attribute</div>
      <div className='options-tray tray'>
        <ShrinkingLabelTextField
          id='attribute-name'
          label='Name (snake_case)'
          value={name}
          onChange={(e: React.ChangeEvent<any>) => setName(e.target.value)}
          pattern={SNAKE_CASE}
        />
        <ShrinkingLabelTextField
          id='attribute-description'
          value={description}
          label='Description (optional)'
          onChange={(e: React.ChangeEvent<any>) =>
            setDescription(e.target.value)
          }
        />
        <ShrinkingLabelTextField
          id='attribute-group'
          value={group}
          label='Group (optional; comma-separated list)'
          onChange={(e: React.ChangeEvent<any>) => setGroup(e.target.value)}
          pattern={COMMA_SEP}
        />
        <ModalSelect
          id='attribute-type'
          value={type}
          label='Type'
          onChange={setType}
          options={attributeTypes}
        />
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
