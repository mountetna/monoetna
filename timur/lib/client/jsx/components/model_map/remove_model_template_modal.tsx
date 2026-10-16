import React from 'react';

import Typography from '@material-ui/core/Typography';

import ModelActionsModal, { ModelModalParams } from './model_actions_modal';

export default function RemoveModelTemplateModal({
  onSave,
  onClose,
  open,
  modelName,
  templateName
}: ModelModalParams & { modelName: string, templateName: string }) {
  return (
    <ModelActionsModal onClose={onClose} open={open} onSave={onSave} title='Remove Template' saveDisabled={false} saveLabel='Remove'>
        <Typography>
          {modelName} will no longer be associated to {templateName}.
        </Typography>
    </ModelActionsModal>
  );
}
