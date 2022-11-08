import React, {useState, useCallback, useEffect, useMemo} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {useModal} from 'etna-js/components/ModalDialogContainer';

import {selectModels} from 'etna-js/selectors/magma';
import {useReduxState} from 'etna-js/hooks/useReduxState';

import {SNAKE_CASE, SNAKE_CASE_STRICT} from '../../utils/edit_map';
import DisabledButton from '../search/disabled_button';
import ModalSelect from './modal_select';
import {ShrinkingLabelTextField} from './shrinking_label_text_field';

export default function AddModelModal({onSave}: {onSave: any}) {
  const [disabled, setDisabled] = useState(true);
  const [identifier, setIdentifier] = useState('');
  const [childModelName, setChildModelName] = useState('');
  const [modelLinkType, setModelLinkType] = useState('');
  const [childModelNameExists, setChildModelNameExists] = useState(false);

  const models = useReduxState((state: any) => selectModels(state));

  const {dismissModal} = useModal();
  const invoke = useActionInvoker();

  const handleOnSave = useCallback(() => {
    onSave({
      identifier,
      model_name: childModelName,
      parent_link_type: modelLinkType
    });
  }, [identifier, childModelName, modelLinkType]);

  const existingModelNames = useMemo(() => {
    return Object.keys(models);
  }, [models]);

  useEffect(() => {
    if (
      identifier &&
      childModelName &&
      modelLinkType &&
      !childModelNameExists
    ) {
      setDisabled(false);
    } else {
      setDisabled(true);
    }
  }, [childModelNameExists, identifier, childModelName, modelLinkType]);

  const handleOnCancel = useCallback(() => {
    invoke(dismissModal());
  }, [invoke, dismissModal]);

  const validateChildModelName = useCallback(
    (input: string) => {
      if (existingModelNames.includes(input)) {
        setChildModelNameExists(true);
      } else {
        setChildModelNameExists(false);
      }
      setChildModelName(input);
    },
    [existingModelNames]
  );

  const modelLinkTypeOptions = ['child', 'collection', 'table'];

  return (
    <div className='add-model-modal model-actions-modal'>
      <div className='header'>Add Model</div>
      <div className='options-tray tray'>
        <ShrinkingLabelTextField
          id='model-identifier'
          label='Identifier (snake_case)'
          value={identifier}
          onChange={(e: React.ChangeEvent<any>) =>
            setIdentifier(e.target.value)
          }
          pattern={SNAKE_CASE}
        />
        <ShrinkingLabelTextField
          id='child-model-name'
          label='Child Model Name (snake_case no numbers)'
          value={childModelName}
          additionalError={childModelNameExists}
          onChange={(e: React.ChangeEvent<any>) =>
            validateChildModelName(e.target.value)
          }
          pattern={SNAKE_CASE_STRICT}
        />
        <ModalSelect
          id='model-type'
          value={modelLinkType}
          label='Link Type'
          onChange={setModelLinkType}
          options={modelLinkTypeOptions}
        />
      </div>
      <div className='options-action-wrapper'>
        <DisabledButton
          id='cancel-add-link-btn'
          className='cancel'
          label='Cancel'
          disabled={false}
          onClick={handleOnCancel}
        />
        <DisabledButton
          id='add-link-btn'
          className='save'
          label='Save'
          disabled={disabled}
          onClick={handleOnSave}
        />
      </div>
    </div>
  );
}
