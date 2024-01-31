import React, {useState, useCallback, useEffect} from 'react';

import {ShrinkingLabelTextField} from './shrinking_label_text_field';
import {COMMA_SEP, SNAKE_CASE} from '../../utils/edit_map';
import ModalSelect from './modal_select';
import ModelActionsModal from './model_actions_modal';

export default function AddAttributeModal({onClose,open,onSave}: {onSave: any}) {
  const [name, setName] = useState('');
  const [description, setDescription] = useState('');
  const [type, setType] = useState('');
  const [group, setGroup] = useState('');

  const handleOnSave = useCallback(() => {
    onSave({
      attribute_name: name,
      description,
      type,
      attribute_group: group
    });
  }, [name, description, type, group]);

  const disabled = !(name && type);

  const handleOnCancel = useCallback(() => {
    onClose();
    reset();
  }, []);

  const reset = useCallback(() => {
    setDisabled(true);
    setName('');
    setDescription('');
    setType('');
    setGroup('');
  }, []);

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
    <ModelActionsModal onClose={handleOnCancel} open={open} onSave={handleOnSave} title='Add Attribute' saveDisabled={disabled}>
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
    </ModelActionsModal>
  );
}
