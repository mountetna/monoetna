import React, {useState, useEffect, useCallback} from 'react';
import Grid from '@material-ui/core/Grid';
import {makeStyles} from '@material-ui/core/styles';
import TextField from '@material-ui/core/TextField';
import InputAdornment from '@material-ui/core/InputAdornment';
import Search from '@material-ui/icons/Search';

import {fetchUsers} from '../api/janus_api';

import {UserFlagsInterface} from '../models/user_models';
import UserTable from './flags-user-table';
import TableControls from './flags-table-controls';

const useStyles = makeStyles((theme) => ({
  margin: {
    margin: theme.spacing(3)
  }
}));

const FlagsView = () => {
  const [allUsers, setAllUsers] = useState([] as UserFlagsInterface[]);
  const [filteredUsers, setFilteredUsers] = useState(
    [] as UserFlagsInterface[]
  );
  const [searchTerm, setSearchTerm] = useState('');
  const [searchProjects, setSearchProjects] = useState([] as string[]);
  const [searchFlags, setSearchFlags] = useState([] as string[]);

  const classes = useStyles();

  useEffect(() => {
    fetchUsers().then(({users}) => {
      setAllUsers(users);
    });
  }, []);

  useEffect(() => {
    setFilteredUsers(allUsers);
  }, [allUsers]);

  useEffect(() => {
    // (searchTerm across user.name || user.email) &&
    //  (searchProjects OR'd) && (searchFlags OR'd)
    setFilteredUsers(
      allUsers
        .filter((user) => {
          let regex = new RegExp(searchTerm);

          return regex.test(user.name) || regex.test(user.email);
        })
        .filter((user) => {
          if (searchProjects.length === 0) return true;

          return user.projects.some((p) => searchProjects.includes(p));
        })
        .filter((user) => {
          if (searchFlags.length === 0) return true;
          if (null == user.flags) return false;

          return user.flags.some((f) => searchFlags.includes(f));
        })
    );
  }, [searchTerm, searchProjects, searchFlags]);

  return (
    <Grid container xs={12} direction='column' className={classes.margin}>
      <Grid item>
        <TableControls
          onChangeSearch={setSearchTerm}
          onChangeProjects={setSearchProjects}
          onChangeFlags={setSearchFlags}
          projectOptions={[...new Set(allUsers.map((u) => u.projects).flat())]}
          flagOptions={[...new Set(allUsers.map((u) => u.flags || []).flat())]}
        />
      </Grid>
      <Grid item>
        <UserTable users={filteredUsers} />
      </Grid>
    </Grid>
  );
};

export default FlagsView;
