import React, {useState, useCallback, useEffect, useMemo} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {useModal} from 'etna-js/components/ModalDialogContainer';

import {selectModels} from 'etna-js/selectors/magma';
import {useReduxState} from 'etna-js/hooks/useReduxState';

import DisabledButton from 'etna-js/components/disabled_button';
import ModalSelect from './modal_select';

export default function ReparentModelModal({
  onSave,
  modelName
}: {
  onSave: any;
  modelName: string;
}) {
  const [disabled, setDisabled] = useState(true);
  const [parentModelName, setParentModelName] = useState('');

  const models = useReduxState((state: any) => selectModels(state));

  const {dismissModal} = useModal();
  const invoke = useActionInvoker();

  const handleOnSave = useCallback(() => {
    onSave(parentModelName);
  }, [parentModelName]);

  const currentParentModelName = useMemo(() => {
    return models[modelName].template.parent;
  }, [modelName, models]);

  const existingModelNames = useMemo(() => {
    return Object.keys(models);
  }, [models]);

  const parentModelNameExists = useMemo(() => {
    return existingModelNames.includes(parentModelName);
  }, [existingModelNames, parentModelName]);

  useEffect(() => {
    if (parentModelName && parentModelNameExists) {
      setDisabled(false);
    } else {
      setDisabled(true);
    }
  }, [parentModelName, parentModelNameExists]);

  const handleOnCancel = useCallback(() => {
    invoke(dismissModal());
  }, [invoke, dismissModal]);

  const parentModelNameOptions = useMemo(() => {
    return existingModelNames.filter((modelName: string) => {
      return modelName !== currentParentModelName;
    });
  }, [existingModelNames, currentParentModelName]);

  return (
    <div className='reparent-model-modal model-actions-modal'>
      <div className='header'>Reparent Model</div>
      <div className='options-tray tray'>
        <ModalSelect
          id='model-type'
          value={parentModelName}
          label='New Parent Model'
          onChange={setParentModelName}
          options={parentModelNameOptions}
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
