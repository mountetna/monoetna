import React, {useState, useCallback, useEffect} from 'react';

import TextField from '@material-ui/core/TextField';
import Checkbox from '@material-ui/core/Checkbox';
import FormControlLabel from '@material-ui/core/FormControlLabel';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {useModal} from 'etna-js/components/ModalDialogContainer';

import DisabledButton from '../search/disabled_button';
import {Attribute} from '../../api/magma_api';

export default function EditAttributeModal({
  onSave,
  attribute
}: {
  onSave: any;
  attribute: Attribute;
}) {
  const [disabled, setDisabled] = useState(true);
  const [updatedAttribute, setUpdatedAttribute] = useState({...attribute});
  const {dismissModal} = useModal();
  const invoke = useActionInvoker();

  const handleOnSave = useCallback(() => {
    onSave({
      ...updatedAttribute,
      attribute_name: attribute.attribute_name,
      new_attribute_name: updatedAttribute.attribute_name
    });
  }, [attribute, updatedAttribute]);

  useEffect(() => {
    if (updatedAttribute.attribute_name) {
      setDisabled(false);
    } else {
      setDisabled(true);
    }
  }, [updatedAttribute]);

  const handleOnCancel = useCallback(() => {
    invoke(dismissModal());
  }, [invoke, dismissModal]);

  const updateAttribute = useCallback(
    (updatePairs: [string, string | boolean][]) => {
      let tmpAttribute: Attribute = {...updatedAttribute};
      updatePairs.forEach(([key, value]) => {
        (tmpAttribute as Record<keyof Attribute, string | boolean>)[
          key as keyof Attribute
        ] = value;
      });
      setUpdatedAttribute(tmpAttribute);
    },
    [updatedAttribute]
  );

  return (
    <div className='edit-attribute-modal model-actions-modal'>
      <div className='header'>Edit Attribute</div>
      <div className='options-tray tray'>
        <TextField
          id='edit-attribute-name'
          label='Name'
          value={updatedAttribute.attribute_name}
          onChange={(e) =>
            updateAttribute([
              ['attribute_name', e.target.value],
              ['new_attribute_name', e.target.value]
            ])
          }
        />
        <TextField
          id='edit-attribute-description'
          label='Description'
          value={updatedAttribute.description}
          onChange={(e) => updateAttribute([['description', e.target.value]])}
        />
        <TextField
          id='edit-attribute-group'
          label='Group (comma-separated list)'
          value={updatedAttribute.attribute_group}
          onChange={(e) =>
            updateAttribute([['attribute_group', e.target.value]])
          }
        />
        <TextField
          id='edit-attribute-display-name'
          label='Display Name'
          value={updatedAttribute.display_name}
          onChange={(e) => updateAttribute([['display_name', e.target.value]])}
        />
        <TextField
          id='edit-attribute-validation'
          label='Validation (JSON object with `type` and `value`)'
          value={JSON.stringify(updatedAttribute.validation)}
          onChange={(e) =>
            updateAttribute([['validation', JSON.parse(e.target.value)]])
          }
        />
        <FormControlLabel
          control={
            <Checkbox
              onChange={(e) => updateAttribute([['hidden', e.target.checked]])}
              checked={updatedAttribute.hidden}
              inputProps={{'aria-label': 'controlled'}}
            />
          }
          label='Hidden'
        />
        <FormControlLabel
          control={
            <Checkbox
              onChange={(e) =>
                updateAttribute([['read_only', e.target.checked]])
              }
              checked={updatedAttribute.read_only}
              inputProps={{'aria-label': 'controlled'}}
            />
          }
          label='Read-only'
        />
        <FormControlLabel
          control={
            <Checkbox
              onChange={(e) =>
                updateAttribute([['restricted', e.target.checked]])
              }
              checked={updatedAttribute.restricted}
              inputProps={{'aria-label': 'controlled'}}
            />
          }
          label='Restricted'
        />
      </div>
      <div className='options-action-wrapper'>
        <DisabledButton
          id='cancel-edit-attribute-btn'
          className='cancel'
          label='Cancel'
          disabled={false}
          onClick={handleOnCancel}
        />
        <DisabledButton
          id='edit-attribute-btn'
          className='save'
          label='Save'
          disabled={disabled}
          onClick={handleOnSave}
        />
      </div>
    </div>
  );
}
