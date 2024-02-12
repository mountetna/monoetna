import React, {useState, useCallback, useEffect, useMemo} from 'react';

import {makeStyles} from '@material-ui/core/styles';
import Typography from '@material-ui/core/Typography';
import DialogContentText from '@material-ui/core/DialogContentText';

import {selectModels} from 'etna-js/selectors/magma';
import {useReduxState} from 'etna-js/hooks/useReduxState';

import {SNAKE_CASE, SNAKE_CASE_STRICT} from '../../utils/edit_map';
import {ShrinkingLabelTextField} from './shrinking_label_text_field';
import AntSwitch from '../query/ant_switch';
import ModelActionsModal, { ModelModalParams } from './model_actions_modal';

const useStyles = makeStyles((theme) => ({
  switch: {
    marginTop: '1rem'
  },
  filler: {
    height: '48px'
  }
}));

export default function AddModelModal({
  modelName,
  onSave,
  onClose,
  open
}: ModelModalParams & { modelName: string }) {
  const [identifier, setIdentifier] = useState('');
  const [childModelName, setChildModelName] = useState('');
  const [isTable, setIsTable] = useState(false);
  const [childModelNameExists, setChildModelNameExists] = useState(false);
  const classes = useStyles();

  const models = useReduxState((state: any) => selectModels(state));

  const handleOnSave = useCallback(() => {
    onSave({
      identifier: isTable ? 'id' : identifier,
      model_name: childModelName,
      parent_link_type: isTable ? 'table' : 'collection'
    });
  }, [identifier, childModelName, isTable]);

  const handleOnCancel = useCallback(() => {
    onClose();
    reset();
  }, []);

  const reset = useCallback(() => {
    setIdentifier('');
    setChildModelName('');
    setIsTable(false);
    setChildModelNameExists(false);
  }, []);

  const existingModelNames = useMemo(() => {
    return Object.keys(models);
  }, [models]);

  const disabled = !((isTable ? true : identifier)
    && childModelName
    && !childModelNameExists);

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

  return (
    <ModelActionsModal onClose={handleOnCancel} open={open} onSave={handleOnSave} title='Add Model' saveDisabled={disabled}>
        <DialogContentText>
          <Typography color='primary'>Add a child model below parent model <Typography component='span' color='secondary'>{modelName}</Typography>.</Typography>
        </DialogContentText>
        <ShrinkingLabelTextField
          id='child-model-name'
          label='New Model Name (snake_case no numbers)'
          value={childModelName}
          additionalError={childModelNameExists}
          onChange={(e: React.ChangeEvent<any>) =>
            validateChildModelName(e.target.value)
          }
          pattern={SNAKE_CASE_STRICT}
        />
        <div className={classes.switch}>
          <AntSwitch
            leftOption='Collection'
            rightOption='Table'
            name='Model Type'
            checked={isTable}
            onChange={() => setIsTable(!isTable)}
          />
          {isTable ? (
            <Typography>
              The model will contain tabular data associated with the parent
              record.
            </Typography>
          ) : (
            <Typography>
              The model will be a collection of named records.
            </Typography>
          )}
        </div>
        {!isTable ? (
          <ShrinkingLabelTextField
            id='model-identifier'
            label='Identifier attribute name (snake_case)'
            value={identifier}
            onChange={(e: React.ChangeEvent<any>) =>
              setIdentifier(e.target.value)
            }
            pattern={SNAKE_CASE}
          />
        ) : (
          <div className={classes.filler}></div>
        )}
    </ModelActionsModal>
  );
}
