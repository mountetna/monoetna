import React from 'react';
import Button from '@material-ui/core/Button';
import TextField from '@material-ui/core/TextField';
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import DialogContentText from '@material-ui/core/DialogContentText';
import DialogTitle from '@material-ui/core/DialogTitle';

export default function ConfirmationCodeModal({open, confirmation, confirmationError, onChange, onClose, onConfirm}:{
  open: boolean,
  confirmation: string,
  confirmationError: string,
  onClose: () => void,
  onChange: (value: string) => void,
  onConfirm: () => void
}) {
  return (
    <Dialog
      open={open}
      onClose={onClose}
      aria-labelledby="message-modal-title"
      aria-describedby="message-modal-description"
    >
      <DialogTitle id="message-modal-title">Confirmation Required</DialogTitle>
      <DialogContent>
        <DialogContentText id="message-modal-description">
          { confirmationError }
        </DialogContentText>
        <TextField
          autoFocus
          required
          margin="dense"
          name="confirmation"
          label="Confirmation"
          type="string"
          fullWidth
          variant="standard"
          onChange={ e => onChange(e.target.value) }
          value={ confirmation }
        />
      </DialogContent>
      <DialogActions>
        <Button onClick={onClose} color="primary" autoFocus>
          Close
        </Button>
        <Button onClick={onConfirm} color="primary" autoFocus>
          Confirm
        </Button>
      </DialogActions>
    </Dialog>
  );
}
