import React from 'react';
import Button from '@material-ui/core/Button';
import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import DialogContentText from '@material-ui/core/DialogContentText';
import DialogTitle from '@material-ui/core/DialogTitle';

export default function MessageModal({open, title, message, onClose}:{ open: boolean, title: string, message: string, onClose: Function }) {
  return (
      <Dialog
        open={open}
        onClose={onClose}
        aria-labelledby="message-modal-title"
        aria-describedby="message-modal-description"
      >
        <DialogTitle id="message-modal-title">{ title }</DialogTitle>
        <DialogContent>
          <DialogContentText id="message-modal-description">
            { message }
          </DialogContentText>
        </DialogContent>
        <DialogActions>
          <Button onClick={onClose} color="primary" autoFocus>
            Close
          </Button>
        </DialogActions>
      </Dialog>
  );
}
