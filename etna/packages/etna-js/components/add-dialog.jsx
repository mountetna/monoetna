import React, {useState, useEffect, useCallback, useReducer} from 'react';

import Dialog from '@material-ui/core/Dialog';
import DialogActions from '@material-ui/core/DialogActions';
import DialogContent from '@material-ui/core/DialogContent';
import DialogContentText from '@material-ui/core/DialogContentText';
import DialogTitle from '@material-ui/core/DialogTitle';
import Button from '@material-ui/core/Button';
import TextField from '@material-ui/core/TextField';
import AddIcon from '@material-ui/icons/Add';

const AddDialog = ({update, buttonClass, title, content, buttonText, placeholders, mask=(e => e), mask2=(e => e)}) => {
  const [ open, setOpen ] = useState(false);
  const [ v1, setV1 ] = useState('');
  const [ v2, setV2 ] = useState('');

  const handleClose = () => {
    setOpen(false);
    setV1('');
    setV2('');
  };

  const handleAdd = useCallback( () => {
    update(v1,v2);
    handleClose();
  }, [v1, v2]);

  return <React.Fragment>
    <Button 
    variant='text'
    startIcon={<AddIcon/>}
    color='secondary'
    className={buttonClass}
    onClick={() => setOpen(true)}>{buttonText}</Button>
    <Dialog open={open} onClose={handleClose} aria-labelledby='form-dialog-title'>
      <DialogTitle id='form-dialog-title'>{title}</DialogTitle>
      <DialogContent>
        <DialogContentText>
          {content}
        </DialogContentText>
        <TextField
          autoFocus
          margin='dense'
          placeholder={ placeholders[0] }
          fullWidth
          value={ v1 }
          onChange={ e => setV1(mask(e.target.value)) }
        />
        { placeholders[1] && <TextField
          margin='dense'
          placeholder={ placeholders[1] }
          fullWidth
          value={ v2 }
          onChange={ e => setV2(mask2(e.target.value)) }
        /> }
      </DialogContent>
      <DialogActions>
        <Button onClick={handleClose} color='primary'>
          Cancel
        </Button>
        <Button onClick={ handleAdd } color='primary' disabled={ !v1 || (placeholders[1] && !v2) }>
          Add
        </Button>
      </DialogActions>
    </Dialog>
  </React.Fragment>;
};


export default AddDialog;
