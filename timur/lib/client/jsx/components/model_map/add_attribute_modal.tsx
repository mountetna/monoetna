import React, {useState, useCallback, useEffect} from 'react';

import {ShrinkingLabelTextField} from './shrinking_label_text_field';
import {COMMA_SEP, SNAKE_CASE, MATRIX_VALIDATION_TYPES} from '../../utils/edit_map';
import ModalSelect from './modal_select';
import ModelActionsModal, { ModelModalParams } from './model_actions_modal';

export default function AddAttributeModal({onClose,open,onSave}: ModelModalParams) {
  const [name, setName] = useState('');
  const [description, setDescription] = useState('');
  const [type, setType] = useState('');
  const [group, setGroup] = useState('');
  const [validationValue, setValidationValue] = useState('');
  const [validationType, setValidationType] = useState('CopyExternal');

  const handleOnSave = useCallback(() => {
    const params: {[k: string]: any} = {
      attribute_name: name,
      description,
      type,
      attribute_group: group
    };

    if (type=='matrix') {
      params['validation'] = {
        type: validationType,
        value: validationValue.split(',').map((s) => s.trim())
      };
    }

    onSave(params);
  }, [name, description, type, group, validationValue, validationType]);

  const disabled = !(name && type) || (type=='matrix' && !validationValue);

  const reset = useCallback(() => {
    setName('');
    setDescription('');
    setType('');
    setGroup('');
    setValidationValue('');
    setValidationType('CopyExternal');
  }, []);

  const handleOnCancel = useCallback(() => {
    onClose();
    reset();
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
    'integer',
    'matrix'
  ];

  function updateValidationType(type) {
    setValidationValue('');
    setValidationType(type);
  }

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
        {type=='matrix' && <ModalSelect
          id='attribute-validation-type'
          value={validationType}
          label='Validation Type'
          onChange={updateValidationType}
          options={MATRIX_VALIDATION_TYPES}
        />}
        {type=='matrix' && <ShrinkingLabelTextField
          id='attribute-validation-value'
          value={validationValue}
          label={validationType === 'Array' ? 
            `Feature Names (comma-separated list) -- UI length limit: ${validationValue.length} / 2000 characters` :
            'Copy Source: <project_name>, <model_name>, <attribute_name>'
          }
          onChange={(e: React.ChangeEvent<any>) => {
            setValidationValue(e.target.value.length <= 2000 ? e.target.value : '')
          }}
          pattern={COMMA_SEP}
        />}
    </ModelActionsModal>
  );
}
