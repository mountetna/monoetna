import React, {useState, useEffect, useCallback, useRef} from 'react';
import Button from '@material-ui/core/Button';
import FormControl from '@material-ui/core/FormControl';
import Alert from '@material-ui/lab/Alert';
import InputLabel from '@material-ui/core/InputLabel';
import Select from '@material-ui/core/Select';
import MenuItem from '@material-ui/core/MenuItem';
import TextField from '@material-ui/core/TextField';
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import DialogContentText from '@material-ui/core/DialogContentText';
import DialogTitle from '@material-ui/core/DialogTitle';
import {json_post} from 'etna-js/utils/fetch';
import {Job, EtnaError} from '../polyphemus';

const WorkflowCreate = ({
  project_name,
  open,
  onClose,
  onCreate,
  jobs
}: {
  project_name: string;
  open: boolean;
  onClose: Function;
  onCreate: Function;
  jobs: Job[];
}) => {
  const [workflow_name, setWorkflowName] = useState('');
  const [workflow_type, setWorkflowType] = useState('');
  const [error, setError] = useState<string | null>(null);

  const reset = () => {
    setWorkflowName('');
    setWorkflowType('');
    setError(null);
  };

  const close = () => {
    reset();
    onClose();
  };

  const createJob = () =>
    json_post(`/api/workflows/${project_name}/create`, {workflow_name, workflow_type})
      .then((workflow) => {
        onCreate(workflow);
        close();
      })
      .catch((response) =>
        response.then(({error}: EtnaError) => setError(error))
      );

  return (
    <Dialog open={open} onClose={close} aria-labelledby='form-dialog-title'>
      <DialogTitle id='form-dialog-title'>Add Loader</DialogTitle>
      <DialogContent>
        {error && <Alert>{error}</Alert>}
        <FormControl fullWidth>
          <TextField
            autoFocus
            margin='dense'
            label='Loader Name'
            value={workflow_name}
            onChange={(e) => setWorkflowName(e.target.value as string)}
          />
        </FormControl>
        <FormControl fullWidth>
          <InputLabel>Job Type</InputLabel>
          <Select
            value={workflow_type}
            onChange={(e) => setWorkflowType(e.target.value as string)}
          >
            {jobs.map(({name}) => (
              <MenuItem key={name} value={name}>
                {name}
              </MenuItem>
            ))}
          </Select>
        </FormControl>
      </DialogContent>
      <DialogActions>
        {workflow_name && workflow_type && (
          <Button onClick={createJob} color='primary'>
            Create
          </Button>
        )}
        <Button onClick={close} color='secondary'>
          Cancel
        </Button>
      </DialogActions>
    </Dialog>
  );
};

export default WorkflowCreate;
