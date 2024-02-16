import React, {useState, useCallback, useEffect, useMemo} from 'react';

import {makeStyles} from '@material-ui/core/styles';
import Typography from '@material-ui/core/Typography';
import Grid from '@material-ui/core/Grid';
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import DialogTitle from '@material-ui/core/DialogTitle';
import Button from '@material-ui/core/Button';

import {selectModels} from 'etna-js/selectors/magma';
import {useReduxState} from 'etna-js/hooks/useReduxState';

import {SNAKE_CASE, SNAKE_CASE_STRICT} from '../../utils/edit_map';
import {ShrinkingLabelTextField} from './shrinking_label_text_field';
import AntSwitch from '../query/ant_switch';
import {SelectProjectModel} from '../select_project_model';

const useStyles = makeStyles((theme) => ({
  switch: {
    marginTop: '1rem'
  },
  filler: {
    height: '48px'
  }
}));

export type ModelModalParams = {
  title: string,
  saveDisabled: boolean;
  saveLabel?: string;
  open: boolean,
  onSave: any,
  onClose: any,
};

export default function ModelActionsModal ({
  onSave,
  title,
  saveDisabled,
  saveLabel='Save',
  open,
  onClose,
  children
}:ModelModalParams & {children: any}) {
  const classes = useStyles();

  const handleSave = () => {
    onSave();
    onClose();
  }

  return (
    <Dialog maxWidth='md' fullWidth open={open}>
      <DialogTitle>{title}</DialogTitle>
      <DialogContent>
        {children}
      </DialogContent>
      <DialogActions>
        <Button color='secondary' onClick={onClose}>Cancel</Button>
        <Button color='primary' disabled={saveDisabled} onClick={handleSave}>{saveLabel}</Button>
      </DialogActions>
    </Dialog>
  );
}
