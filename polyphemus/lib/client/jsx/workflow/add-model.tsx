import React, { useState, useContext } from 'react';
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import DialogContentText from '@material-ui/core/DialogContentText';
import DialogTitle from '@material-ui/core/DialogTitle';
import Select from '@material-ui/core/Select';
import MenuItem from '@material-ui/core/MenuItem';
import Button from '@material-ui/core/Button';
import {MagmaContext} from 'etna-js/contexts/magma-context';
import {diff} from '../utils/list';

const AddModel = ({
  open,
  close,
  update,
  excludeModels
}: {
  open: boolean;
  close: () => void;
  update: (newModelName: string) => void;
  excludeModels: string[];
}) => {
  const [newModel, setNewModel] = useState('');

  const {models} = useContext(MagmaContext);

  const model_names = diff(Object.keys(models), excludeModels).sort();

  return (
    <Dialog open={open} onClose={close}>
      <DialogTitle>Add Model</DialogTitle>
      <DialogContent>
        <Select
          displayEmpty
          value={newModel}
          onChange={(e) => setNewModel(e.target.value as string)}
        >
          <MenuItem value=''>
            <em>None</em>
          </MenuItem>
          {model_names.map((att_name) => (
            <MenuItem key={att_name} value={att_name}>
              {att_name}
            </MenuItem>
          ))}
        </Select>
      </DialogContent>
      <DialogActions>
        <Button
          disabled={!newModel}
          onClick={() => {
            update(newModel);
            close();
          }}
          color='secondary'
        >
          Add
        </Button>
      </DialogActions>
    </Dialog>
  );
};

export default AddModel
