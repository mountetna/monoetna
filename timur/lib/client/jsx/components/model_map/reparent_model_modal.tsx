import React, {useState, useCallback, useEffect, useMemo} from 'react';

import {selectModels} from 'etna-js/selectors/magma';
import {useReduxState} from 'etna-js/hooks/useReduxState';

import ModalSelect from './modal_select';
import ModelActionsModal, { ModelModalParams } from './model_actions_modal';

export default function ReparentModelModal({
  onSave,
  onClose,
  open,
  modelName
}: ModelModalParams & { modelName: string }) {
  const [parentModelName, setParentModelName] = useState('');

  const models = useReduxState((state: any) => selectModels(state));

  const handleOnSave = useCallback(() => {
    onSave(parentModelName);
  }, [parentModelName]);

  const currentParentModelName = useMemo(() => {
    return models[modelName]?.template?.parent;
  }, [modelName, models]);

  const existingModelNames = useMemo(() => {
    return Object.keys(models);
  }, [models]);

  const disabled = !(parentModelName
    && existingModelNames.includes(parentModelName));

  const handleOnCancel = useCallback(() => {
    onClose();
    reset();
  }, []);

  const reset = useCallback(() => {
    setParentModelName('');
  }, []);

  const parentModelNameOptions = useMemo(() => {
    return existingModelNames.filter((modelName: string) => {
      return modelName !== currentParentModelName;
    });
  }, [existingModelNames, currentParentModelName]);

  return (
    <ModelActionsModal onClose={handleOnCancel} open={open} onSave={handleOnSave} title='Reparent Model' saveDisabled={disabled}>
        <ModalSelect
          id='model-type'
          value={parentModelName}
          label='New Parent Model'
          onChange={setParentModelName}
          options={parentModelNameOptions}
        />
    </ModelActionsModal>
  );
}
