import React, {useState, useCallback, useEffect, useMemo} from 'react';

import {makeStyles} from '@material-ui/core/styles';
import Typography from '@material-ui/core/Typography';

import {useActionInvoker} from 'etna-js/hooks/useActionInvoker';
import {useModal} from 'etna-js/components/ModalDialogContainer';
import {selectModels} from 'etna-js/selectors/magma';
import {useReduxState} from 'etna-js/hooks/useReduxState';

import {SNAKE_CASE, SNAKE_CASE_STRICT} from '../../utils/edit_map';
import DisabledButton from '../search/disabled_button';
import {ShrinkingLabelTextField} from './shrinking_label_text_field';
import AntSwitch from '../query/ant_switch';

const useStyles = makeStyles((theme) => ({
  switch: {
    marginTop: '1rem'
  },
  filler: {
    height: '48px'
  }
}));

export default function AddModelModal({onSave}: {onSave: any}) {
  const [disabled, setDisabled] = useState(true);
  const [identifier, setIdentifier] = useState('');
  const [childModelName, setChildModelName] = useState('');
  const [isTable, setIsTable] = useState(false);
  const [childModelNameExists, setChildModelNameExists] = useState(false);
  const classes = useStyles();

  const models = useReduxState((state: any) => selectModels(state));

  const {dismissModal} = useModal();
  const invoke = useActionInvoker();

  const handleOnSave = useCallback(() => {
    onSave({
      identifier: isTable ? 'id' : identifier,
      model_name: childModelName,
      parent_link_type: isTable ? 'table' : 'collection'
    });
  }, [identifier, childModelName, isTable]);

  const existingModelNames = useMemo(() => {
    return Object.keys(models);
  }, [models]);

  useEffect(() => {
    if (
      (isTable ? true : identifier) &&
      childModelName &&
      !childModelNameExists
    ) {
      setDisabled(false);
    } else {
      setDisabled(true);
    }
  }, [childModelNameExists, identifier, childModelName, isTable]);

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

  return (
    <div className='add-model-modal model-actions-modal'>
      <div className='header'>Add Model</div>
      <div className='options-tray tray'>
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
