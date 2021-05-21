import React, {useState, useEffect, useContext} from 'react';
import Button from '@material-ui/core/Button';
import {makeStyles} from '@material-ui/core/styles';
import Card from '@material-ui/core/Card';
import CardHeader from '@material-ui/core/CardHeader';
import TextField from '@material-ui/core/TextField';
import CardContent from '@material-ui/core/CardContent';
import CardActions from '@material-ui/core/CardActions';
import AddIcon from '@material-ui/icons/Add';
import RemoveIcon from '@material-ui/icons/Remove';

import {FlagsContext} from './flags-context';
import {UserFlagsInterface} from '../types/janus_types';
import {updateUserFlags, fetchUsers} from '../api/janus_api';

const useStyles = makeStyles((theme) => ({
  button: {
    margin: theme.spacing(1)
  },
  card: {
    position: 'fixed',
    zIndex: 10
  },
  error: {
    color: 'red',
    fontWeight: 'bolder'
  }
}));

interface UpdatePayload {
  email: string;
  flags: string[];
}

const AddRemoveFlag = ({
  selectedUsers
}: {
  selectedUsers: UserFlagsInterface[];
}) => {
  const [inputFlag, setInputFlag] = useState('' as string);
  const [payloads, setPayloads] = useState([] as UpdatePayload[]);
  const [message, setMessage] = useState('' as string);
  const [hasError, setHasError] = useState(false);
  let {setUsers} = useContext(FlagsContext);
  const classes = useStyles();

  useEffect(() => {
    if (payloads.length > 0) {
      executeUpdate(payloads);
    }
  }, [payloads]);

  function reset() {
    setHasError(false);
    setMessage('');
    setPayloads([]);
    setInputFlag('');
  }

  function onAddFlag() {
    if ('' === inputFlag) return;

    setPayloads(
      selectedUsers.map((user) => ({
        email: user.email,
        flags: [...(user.flags || [])].concat([inputFlag])
      }))
    );
  }

  function onRemoveFlag() {
    if ('' === inputFlag) return;

    setPayloads(
      selectedUsers.map((user) => ({
        email: user.email,
        flags: [...(user.flags || [])].filter((f) => f !== inputFlag)
      }))
    );
  }

  function executeUpdate(payloads: UpdatePayload[]) {
    let promises: Promise<any>[] = payloads.map((payload) =>
      updateUserFlags(payload)
    );

    Promise.all(promises)
      .then(() => {
        reset();
        return fetchUsers();
      })
      .then(({users}) => {
        setMessage('Saved!');
        setUsers(users);
      })
      .catch((error) => {
        error.then((err: {error: string}) => {
          setMessage(err.error);
          setHasError(true);
        });
      });
  }

  return (
    <Card className={classes.card}>
      <CardHeader
        title='Add or remove a flag'
        subheader={`on ${selectedUsers.length} user${
          selectedUsers.length > 1 ? 's' : ''
        }`}
      />
      <CardContent>
        <TextField
          error={hasError}
          label='Flag text'
          variant='outlined'
          value={inputFlag}
          helperText={message ? message : ''}
          onChange={(e) => setInputFlag(e.target.value as string)}
        />
      </CardContent>
      <CardActions disableSpacing>
        <Button
          variant='contained'
          color='primary'
          className={classes.button}
          startIcon={<AddIcon />}
          onClick={onAddFlag}
        >
          Add Flag
        </Button>
        <Button
          variant='contained'
          color='secondary'
          className={classes.button}
          startIcon={<RemoveIcon />}
          onClick={onRemoveFlag}
        >
          Remove Flag
        </Button>
      </CardActions>
    </Card>
  );
};

export default AddRemoveFlag;
