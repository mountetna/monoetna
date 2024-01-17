import React, {useState, useCallback, useEffect, useMemo} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {useModal} from 'etna-js/components/ModalDialogContainer';

import {selectModels} from 'etna-js/selectors/magma';
import {useReduxState} from 'etna-js/hooks/useReduxState';

import {SNAKE_CASE} from '../../utils/edit_map';
import ModalSelect from './modal_select';
import ModelActionsModal from './model_actions_modal';
import {ShrinkingLabelTextField} from './shrinking_label_text_field';

export default function AddLinkModal({onSave,open,onClose}: {onSave: any}) {
  const [disabled, setDisabled] = useState(true);
  const [linkAttributeName, setLinkAttributeName] = useState('');
  const [reciprocalModelName, setReciprocalModelName] = useState('');
  const [reciprocalAttributeName, setReciprocalAttributeName] = useState('');
  const [reciprocalLinkType, setReciprocalLinkType] = useState('');

  const models = useReduxState((state: any) => selectModels(state));

  const handleOnSave = useCallback(() => {
    onSave({
      linkAttributeName,
      reciprocalModelName,
      reciprocalAttributeName,
      reciprocalLinkType
    });
  }, [
    linkAttributeName,
    reciprocalModelName,
    reciprocalAttributeName,
    reciprocalLinkType
  ]);

  useEffect(() => {
    if (
      linkAttributeName &&
      reciprocalModelName &&
      reciprocalAttributeName &&
      reciprocalLinkType
    ) {
      setDisabled(false);
    } else {
      setDisabled(true);
    }
  }, [
    linkAttributeName,
    reciprocalModelName,
    reciprocalAttributeName,
    reciprocalLinkType
  ]);

  const handleOnCancel = useCallback(() => {
    onClose()
    setDisabled(true);
    setLinkAttributeName('');
    setReciprocalModelName('');
    setReciprocalAttributeName('');
    setReciprocalLinkType('');
  }, []);

  const reciprocalModelNameOptions = useMemo(() => {
    return Object.keys(models);
  }, [models]);

  const reciprocalLinkTypeOptions = ['child', 'collection'];

  return (
    <ModelActionsModal onClose={handleOnCancel} open={open} onSave={handleOnSave} title='Add Link' saveDisabled={disabled}>
        <ShrinkingLabelTextField
          id='link-attribute-name'
          label='Link Attribute Name (snake_case)'
          value={linkAttributeName}
          onChange={(e: React.ChangeEvent<any>) =>
            setLinkAttributeName(e.target.value)
          }
          pattern={SNAKE_CASE}
        />
        <ModalSelect
          id='reciprocal-model-name'
          value={reciprocalModelName}
          label='Reciprocal Model Name'
          onChange={setReciprocalModelName}
          options={reciprocalModelNameOptions}
        />
        {reciprocalModelName && (
          <ShrinkingLabelTextField
            id='reciprocal-attribute-name'
            label='Reciprocal Attribute Name (snake_case)'
            value={reciprocalAttributeName}
            onChange={(e: React.ChangeEvent<any>) =>
              setReciprocalAttributeName(e.target.value)
            }
            pattern={SNAKE_CASE}
          />
        )}
        <ModalSelect
          id='reciprocal-link-type'
          value={reciprocalLinkType}
          label='Reciprocal Link Type'
          onChange={setReciprocalLinkType}
          options={reciprocalLinkTypeOptions}
        />
    </ModelActionsModal>
  );
}
