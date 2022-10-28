import React, {useState, useCallback, useEffect} from 'react';

import TextField from '@material-ui/core/TextField';
import Checkbox from '@material-ui/core/Checkbox';
import MenuItem from '@material-ui/core/MenuItem';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import {makeStyles} from '@material-ui/core/styles';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {useModal} from 'etna-js/components/ModalDialogContainer';

import DisabledButton from '../search/disabled_button';
import {Attribute} from '../../api/magma_api';
import {SNAKE_CASE, COMMA_SEP, VALIDATION_TYPES} from '../../utils/edit_map';
import {ShrinkingLabelTextField} from './shrinking_label_text_field';

const useStyles = makeStyles((theme) => ({
  popover: {
    zIndex: '30000 !important' as any // etna modal is 20000
  }
}));

export default function EditAttributeModal({
  onSave,
  attribute
}: {
  onSave: any;
  attribute: Attribute;
}) {
  const [disabled, setDisabled] = useState(true);
  const [updatedAttribute, setUpdatedAttribute] = useState({...attribute});
  const [validationType, setValidationType] = useState(
    attribute.validation ? attribute.validation.type : ''
  );
  const [validationValue, setValidationValue] = useState(
    attribute.validation ? attribute.validation.value : ''
  );
  const {dismissModal} = useModal();
  const invoke = useActionInvoker();
  const classes = useStyles();

  const isArrayValidation = 'Array' === validationType;

  const handleOnSave = useCallback(() => {
    let params = {
      ...updatedAttribute,
      attribute_name: attribute.attribute_name,
      new_attribute_name: updatedAttribute.attribute_name
    };

    if (validationType && validationValue) {
      params.validation = {
        type: validationType,
        value: isArrayValidation ? validationValue.split(',') : validationValue
      };
    } else {
      params.validation = null;
    }

    onSave(params);
  }, [
    attribute,
    updatedAttribute,
    validationType,
    validationValue,
    isArrayValidation
  ]);

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
        (
          tmpAttribute as Record<
            keyof Attribute,
            string | boolean | {[key: string]: any}
          >
        )[key as keyof Attribute] = value;
      });

      setUpdatedAttribute(tmpAttribute);
    },
    [updatedAttribute]
  );

  return (
    <div className='edit-attribute-modal model-actions-modal'>
      <div className='header'>Edit Attribute</div>
      <div className='options-tray tray'>
        <ShrinkingLabelTextField
          id='edit-attribute-name'
          label='Name (snake_case)'
          value={updatedAttribute.attribute_name}
          onChange={(e: React.ChangeEvent<any>) =>
            updateAttribute([
              ['attribute_name', e.target.value],
              ['new_attribute_name', e.target.value]
            ])
          }
          pattern={SNAKE_CASE}
        />
        <ShrinkingLabelTextField
          id='edit-attribute-description'
          label='Description'
          value={updatedAttribute.description}
          onChange={(e: React.ChangeEvent<any>) =>
            updateAttribute([['description', e.target.value]])
          }
        />
        <ShrinkingLabelTextField
          id='edit-attribute-group'
          label='Group (comma-separated list)'
          value={updatedAttribute.attribute_group}
          onChange={(e: React.ChangeEvent<any>) =>
            updateAttribute([['attribute_group', e.target.value]])
          }
          pattern={COMMA_SEP}
        />
        <ShrinkingLabelTextField
          id='edit-attribute-display-name'
          label='Display Name'
          value={updatedAttribute.display_name}
          onChange={(e: React.ChangeEvent<any>) =>
            updateAttribute([['display_name', e.target.value]])
          }
        />
        <TextField
          id='edit-attribute-validation-type'
          select
          value={validationType}
          label='Validation Type'
          SelectProps={{
            MenuProps: {
              PopoverClasses: {
                root: classes.popover
              }
            }
          }}
          onChange={(e: any) => {
            setValidationValue('');
            setValidationType(e.target.value);
          }}
        >
          {VALIDATION_TYPES.sort().map((option, i) => (
            <MenuItem key={i} value={option}>
              {option}
            </MenuItem>
          ))}
        </TextField>
        {validationType && (
          <ShrinkingLabelTextField
            id='edit-attribute-validation-value'
            label={`Validation ${
              isArrayValidation ? 'Array (comma-separated list)' : 'Regex'
            }`}
            value={validationValue}
            onChange={(e: React.ChangeEvent<any>) =>
              setValidationValue(e.target.value)
            }
            pattern={isArrayValidation ? COMMA_SEP : null}
          />
        )}
        <FormControlLabel
          control={
            <Checkbox
              onChange={(e: React.ChangeEvent<any>) =>
                updateAttribute([['hidden', e.target.checked]])
              }
              checked={updatedAttribute.hidden}
              inputProps={{'aria-label': 'controlled'}}
            />
          }
          label='Hidden'
        />
        <FormControlLabel
          control={
            <Checkbox
              onChange={(e: React.ChangeEvent<any>) =>
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
              onChange={(e: React.ChangeEvent<any>) =>
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
