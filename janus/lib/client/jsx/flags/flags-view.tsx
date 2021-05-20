import React, {useState, useEffect, useCallback} from 'react';
import ScopedCssBaseline from '@material-ui/core/ScopedCssBaseline';
import {makeStyles} from '@material-ui/core/styles';

import {fetchUsers} from '../api/janus_api';

import {UserFlagsInterface} from '../types/janus_types';
import UserTable from './flags-user-table';

const useStyles = makeStyles((theme) => ({
  margin: {
    margin: theme.spacing(3)
  }
}));

const FlagsView = () => {
  // Really, this should go into a Context.
  //   And then we wouldn't have to pass users down through
  //   props, and trigger data-fetching up via props.
  const [allUsers, setAllUsers] = useState([] as UserFlagsInterface[]);

  const classes = useStyles();

  useEffect(() => {
    fetchUsers().then(({users}) => {
      setAllUsers(users);
    });
  }, []);

  return (
    <ScopedCssBaseline>
      <UserTable
        users={allUsers}
        onChange={() => {
          fetchUsers().then(({users}) => {
            setAllUsers(users);
          });
        }}
      />
    </ScopedCssBaseline>
  );
};

export default FlagsView;
