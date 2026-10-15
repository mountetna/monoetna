import React, {useState, useCallback, useEffect, useMemo} from 'react';

import {getModels} from 'etna-js/api/magma_api';
import {selectModels} from 'etna-js/selectors/magma';
import {useReduxState} from 'etna-js/hooks/useReduxState';

import ModalSelect from './modal_select';
import ModelActionsModal, { ModelModalParams } from './model_actions_modal';

export default function SetModelTemplateModal({
  onSave,
  onClose,
  open,
  modelName
}: ModelModalParams & { modelName: string }) {
  const models = useReduxState((state: any) => selectModels(state));

  const currentTemplateModelName = useMemo(() => {
    return models[modelName]?.template?.template_model_name || '';
  }, [modelName, models]);

  const [templateModelName, setTemplateModelName] = useState(currentTemplateModelName);
  const [templateModelNames, setTemplateModelNames] = useState<string[]>([]);

  useEffect(() => {
    if (!open) return;

    getModels(CONFIG.template_project)
      .then(({models}: {models: {[name: string]: any}}) =>
        setTemplateModelNames(Object.keys(models).filter((name) => name !== 'project'))
      )
      .catch(() => setTemplateModelNames([]));
  }, [open]);

  const reset = useCallback(() => {
    setTemplateModelName(currentTemplateModelName);
  }, [currentTemplateModelName]);

  useEffect(() => {
    reset();
  }, [currentTemplateModelName]);

  const handleOnSave = useCallback(() => {
    onSave(templateModelName);
  }, [templateModelName, onSave]);

  const handleOnCancel = useCallback(() => {
    onClose();
    reset();
  }, [onClose, reset]);

  const disabled = !templateModelName || templateModelName === currentTemplateModelName;

  return (
    <ModelActionsModal onClose={handleOnCancel} open={open} onSave={handleOnSave} title='Map to a template model' saveDisabled={disabled}>
        <ModalSelect
          id='template-model'
          value={templateModelName}
          label={`Template Model (${CONFIG.template_project})`}
          onChange={setTemplateModelName}
          options={templateModelNames}
        />
    </ModelActionsModal>
  );
}
