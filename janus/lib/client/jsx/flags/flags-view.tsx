import React, {useState, useEffect, useCallback} from 'react';
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
  let [allUsers, setAllUsers] = useState([] as UserFlagsInterface[]);
  let [filteredUsers, setFilteredUsers] = useState([] as UserFlagsInterface[]);
  const classes = useStyles();

  useEffect(() => {
    fetchUsers().then(({users}) => {
      setAllUsers(users);
    });
  }, []);

  useEffect(() => {
    setFilteredUsers(allUsers);
  }, [allUsers]);

  return (
    <Grid container xs={12} direction='column' className={classes.margin}>
      <Grid item>
        <TextField
          label='Search'
          variant='outlined'
          InputLabelProps={{
            shrink: true
          }}
          InputProps={{
            startAdornment: (
              <InputAdornment position='start'>
                <Search />
              </InputAdornment>
            )
          }}
        />
      </Grid>
      <Grid item>
        <UserTable users={filteredUsers} />
      </Grid>
    </Grid>
  );
};

export default FlagsView;
