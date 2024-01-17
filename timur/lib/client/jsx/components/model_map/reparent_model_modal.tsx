import React, {useState, useCallback, useEffect, useMemo} from 'react';

import {selectModels} from 'etna-js/selectors/magma';
import {useReduxState} from 'etna-js/hooks/useReduxState';

import ModalSelect from './modal_select';
import ModelActionsModal from './model_actions_modal';

export default function ReparentModelModal({
  onSave,
  onClose,
  open,
  modelName
}: {
  onSave: any;
  modelName: string;
}) {
  const [disabled, setDisabled] = useState(true);
  const [parentModelName, setParentModelName] = useState('');

  const models = useReduxState((state: any) => selectModels(state));

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
    onClose();
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
