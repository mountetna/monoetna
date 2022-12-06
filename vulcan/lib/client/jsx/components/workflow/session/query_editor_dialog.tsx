import React, {useState} from 'react';

import DialogTitle from '@material-ui/core/DialogTitle';
import DialogContent from '@material-ui/core/DialogContent';
import DialogActions from '@material-ui/core/DialogActions';
import {makeStyles} from '@material-ui/core/styles';
import SaveIcon from '@material-ui/icons/Save';
import CancelIcon from '@material-ui/icons/Cancel';
import Button from '@material-ui/core/Button';

import QueryPage from 'etna-js/components/query/query_page';
import VulcanQueryControlButtons from 'etna-js/components/query/vulcan_query_control_buttons';

const useStyles = makeStyles((theme) => ({
  dialog: {
    maxWidth: '90vw',
    maxHeight: '90vh'
  },
  subtitle: {display: 'inline'},
  button: {
    margin: '1rem'
  },
  loading: {
    marginLeft: '1rem'
  },
  helpdoc: {
    maxWidth: '600px',
    marginTop: '1rem',
    marginBottom: '1rem'
  },
  propagateButton: {
    marginBottom: '1rem'
  }
}));

function QueryEditorDialog({
  onClose,
  onSave
}: {
  onClose: () => void;
  onSave: ({
    query,
    userColumns
  }: {
    query: string | any[];
    userColumns: string[];
  }) => void;
}) {
  const [query, setQuery] = useState('');
  const [userColumns, setUserColumns] = useState([]);
  const classes = useStyles();

  return (
    <>
      <DialogTitle>Edit Query</DialogTitle>
      <DialogContent className={classes.dialog}>
        <QueryPage
          syncQueryParams={false}
          queryControlButtons={
            <VulcanQueryControlButtons
              onUpdateQuery={({query, userColumns}) => {
                setQuery(query);
                setUserColumns(userColumns);
              }}
            />
          }
        />
      </DialogContent>
      <DialogActions>
        <Button
          onClick={onClose}
          startIcon={<CancelIcon />}
          color='secondary'
          variant='contained'
        >
          Cancel
        </Button>
        <Button
          onClick={() => {
            onSave({
              query,
              userColumns
            });
            onClose();
          }}
          startIcon={<SaveIcon />}
          color='primary'
          variant='contained'
        >
          Save
        </Button>
      </DialogActions>
    </>
  );
}

export default QueryEditorDialog;
