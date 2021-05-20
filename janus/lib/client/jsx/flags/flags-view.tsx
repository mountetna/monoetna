import React, {useState, useEffect, useCallback} from 'react';
import ScopedCssBaseline from '@material-ui/core/ScopedCssBaseline';
import Grid from '@material-ui/core/Grid';
import {makeStyles} from '@material-ui/core/styles';
import TextField from '@material-ui/core/TextField';
import InputAdornment from '@material-ui/core/InputAdornment';
import Search from '@material-ui/icons/Search';

import {fetchUsers} from '../api/janus_api';

import {UserFlagsInterface} from '../models/user_models';
import UserTable from './flags-user-table';

const useStyles = makeStyles((theme) => ({
  margin: {
    margin: theme.spacing(3)
  }
}));

const FlagsView = () => {
  const [allUsers, setAllUsers] = useState([] as UserFlagsInterface[]);

  const classes = useStyles();

  useEffect(() => {
    fetchUsers().then(({users}) => {
      setAllUsers(users);
    });
  }, []);

  return (
    <ScopedCssBaseline>
      <UserTable users={allUsers} />
    </ScopedCssBaseline>
  );
};

export default FlagsView;
