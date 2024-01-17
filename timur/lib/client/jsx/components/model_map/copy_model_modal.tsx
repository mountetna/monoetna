import React, {useState, useCallback, useEffect, useMemo} from 'react';

import {makeStyles} from '@material-ui/core/styles';
import DialogContentText from '@material-ui/core/DialogContentText';
import Typography from '@material-ui/core/Typography';

import {selectModels} from 'etna-js/selectors/magma';
import {useReduxState} from 'etna-js/hooks/useReduxState';

import {SNAKE_CASE, SNAKE_CASE_STRICT} from '../../utils/edit_map';
import {ShrinkingLabelTextField} from './shrinking_label_text_field';
import AntSwitch from '../query/ant_switch';
import {SelectProjectModel} from '../select_project_model';
import ModelActionsModal from './model_actions_modal';

import {useDispatch} from 'react-redux';
import {requestAnswer} from 'etna-js/actions/magma_actions';

const useStyles = makeStyles((theme) => ({
  switch: {
    marginTop: '1rem'
  },
  filler: {
    height: '48px'
  }
}));

export default function CopyModelModal({
  modelName,
  open,
  onSave,
  onClose
}: {
  onSave: any;
  modelName: string;
}) {
  const [disabled, setDisabled] = useState(true);
  const [templateProjectName, setTemplateProjectName] = useState(null);
  const [templateModelName, setTemplateModelName] = useState('');
  const [error, setError] = useState(null);
  const classes = useStyles();

  const models = useReduxState((state: any) => selectModels(state));

  const model = models[templateModelName];

  useEffect( () => {
    if (templateModelName && templateProjectName) {
    }
    requestDocuments({
      project_name: templateProjectName
    })(dispatch)
      .then(({answer}) => {
        setModels(answer.sort());
        setError(null);
      })
      .catch((e) => e.then(({error}) => setError(error)));
  }, [project_name]);

  }, [ templateModelName ]);

  const handleOnSave = useCallback(() => {
    onSave({ });
  }, []);

  const handleOnCancel = useCallback(() => {
    onClose()
    setDisabled(true);
    setTemplateProjectName(null);
    setTemplateModelName('');
    setError(null);
  }, []);

  const dispatch = useDispatch();

  return (
    <ModelActionsModal onClose={handleOnCancel} open={open} onSave={handleOnSave} title='Copy Model' saveDisabled={disabled}>
      <SelectProjectModel
        project_name={templateProjectName}
        setProjectName={setTemplateProjectName}
        model_name={templateModelName}
        setModelName={setTemplateModelName}
        error={error}
        setError={setError}
      />
      { model && <>
        <DialogContentText>
          <Typography>Select attributes to copy</Typography>
        </DialogContentText>
      </>}
      
    </ModelActionsModal>
  );
}
