import React, {useState, useCallback, useEffect, useMemo} from 'react';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {useModal} from 'etna-js/components/ModalDialogContainer';

import {makeStyles} from '@material-ui/core/styles';
import Typography from '@material-ui/core/Typography';

import DisabledButton from '../search/disabled_button';
import {ShrinkingLabelTextField} from './shrinking_label_text_field';

const useStyles = makeStyles((theme) => ({
  instructions: {
    marginBottom: '1rem'
  }
}));

export default function RemoveModelModal({
  onSave,
  onClose,
  open,
  modelName
}: {
  onSave: any;
  modelName: string;
}) {
  const [disabled, setDisabled] = useState(true);
  const [deleteModelName, setDeleteModelName] = useState('');

  const classes = useStyles();

  const handleOnSave = useCallback(() => {
    onSave();
  }, []);

  const deleteModelNameMatches = useMemo(() => {
    return modelName === deleteModelName;
  }, [modelName, deleteModelName]);

  useEffect(() => {
    if (!deleteModelNameMatches) {
      setDisabled(true);
    } else {
      setDisabled(false);
    }
  }, [deleteModelNameMatches]);

  const handleOnCancel = useCallback(() => {
    onClose();
    setDeleteModelName('');
    setDisabled(true);
  }, []);

  return (
    <ModelActionsModal onClose={handleOnCancel} open={open} onSave={handleOnSave} title='Remove Model' saveDisabled={disabled} saveLabel='Remove'>
        <Typography className={classes.instructions}>
          Removing the model may result in loss of data -- please type in the
          name of the model to confirm this action. Note that you cannot undo
          this!
        </Typography>
        <ShrinkingLabelTextField
          id='current-model-name'
          label='Selected Model'
          value={modelName}
          disabled={true}
        />
        <ShrinkingLabelTextField
          id='remove-model-name'
          label='Re-type Model Name'
          value={deleteModelName}
          onChange={(e: React.ChangeEvent<any>) =>
            setDeleteModelName(e.target.value)
          }
        />
    </ModelActionsModal>
  );
}
