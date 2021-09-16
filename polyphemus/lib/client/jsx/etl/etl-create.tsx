import React, { useState, useEffect, useCallback, useRef } from 'react';
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
import { json_post } from 'etna-js/utils/fetch';

const EtlCreate = ({project_name, open, onClose, onCreate, jobs}) => {
  const [ job_name, setJobName ] = useState('');
  const [ job_type, setJobType ] = useState('');
  const [ error, setError ] = useState(null);

  const reset = () => {
    setJobName('');
    setJobType('');
    setError(null);
  };

  const close = () => {
    reset();
    onClose();
  };

  const createJob = () => json_post(
    `/api/etl/${project_name}/create/${job_name}`, {job_type}
  ).then(
    etl => {
      onCreate(etl);
      close();
    }
  ).catch( response => response.then( ({error}) => setError(error)))

  return <Dialog open={open} onClose={ close } aria-labelledby="form-dialog-title">
    <DialogTitle id="form-dialog-title">Add Loader</DialogTitle>
    <DialogContent>
      { error && <Alert>{error}</Alert>}
      <FormControl fullWidth>
        <TextField
          autoFocus
          margin="dense"
          label="Loader Name"
          value={ job_name }
          onChange={ e => setJobName(e.target.value) }
          />
      </FormControl>
      <FormControl fullWidth>
        <InputLabel>Job Type</InputLabel>
        <Select
          value={job_type}
          onChange={ e => setJobType(e.target.value) }
        >
          {
            jobs.map( ({name}) => <MenuItem key={name} value={name}>{name}</MenuItem> )
          }
        </Select>
      </FormControl>
    </DialogContent>
    <DialogActions>
      {
        job_name && job_type && <Button onClick={ createJob } color="primary">Create</Button>
      }
      <Button onClick={ close } color="secondary">Cancel</Button>
    </DialogActions>
  </Dialog>;
}

export default EtlCreate;
