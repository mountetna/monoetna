import React, {useState, useEffect, useContext} from 'react';

import {makeStyles} from '@material-ui/core/styles';
import Checkbox from '@material-ui/core/Checkbox';
import Button from '@material-ui/core/Button';
import Grid from '@material-ui/core/Grid';
import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableContainer from '@material-ui/core/TableContainer';
import TableHead from '@material-ui/core/TableHead';
import TableRow from '@material-ui/core/TableRow';
import Paper from '@material-ui/core/Paper';

import {fetchUsers} from '../api/janus_api';
import {FlagsContext} from './flags-context';
import {UserFlagsInterface} from '../types/janus_types';
import UserRow from './flags-user-row';
import TableControls from './flags-table-controls';
import AddRemove from './flags-add-remove';

const useStyles = makeStyles((theme) => ({
  controls: {
    padding: "15px"
  },
  manage: {
    position: 'relative'
  },
  table: {
    minWidth: 650
  }
}));

const UserTable = () => {
  const [filteredUsers, setFilteredUsers] = useState(
    [] as UserFlagsInterface[]
  );
  const [searchTerm, setSearchTerm] = useState('');
  const [searchProjects, setSearchProjects] = useState([] as string[]);
  const [searchFlags, setSearchFlags] = useState([] as string[]);
  const [selected, setSelected] = useState([] as UserFlagsInterface[]);
  const [addRemoveOpen, setAddRemoveOpen] = useState(false);
  const classes = useStyles();

  let {
    state: {users},
    setUsers
  } = useContext(FlagsContext);

  useEffect(() => {
    fetchUsers().then(({users}) => setUsers(users));
  }, []);

  useEffect(() => {
    if (users) setFilteredUsers(users);
  }, [users]);

  useEffect(() => {
    // (searchTerm across user.name || user.email) &&
    //  (searchProjects OR'd) && (searchFlags OR'd)
    if (users) {
      setFilteredUsers(
        users
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
      setSelected([]);
    }
  }, [searchTerm, searchProjects, searchFlags]);

  if (!users) return null;

  function onSelectAllClick() {
    setSelected(selected.length > 0 ? [] : filteredUsers);
  }

  function onClickUser(user: UserFlagsInterface) {
    if (!isSelected(user)) {
      setSelected([...selected].concat([user]));
    } else {
      setSelected([...selected].filter((u) => u.email !== user.email));
    }
  }

  function isSelected(user: UserFlagsInterface) {
    return selected.filter((u) => u.email === user.email).length > 0;
  }

  function openAddRemoveFlags() {
    setAddRemoveOpen(!addRemoveOpen);
  }

  return (
    <Grid container direction='column'>
      <Grid container item xs={12} className={classes.controls}>
        <Grid container item xs={9} alignItems='flex-end' justify='flex-start'>
          <TableControls
            onChangeSearch={setSearchTerm}
            onChangeProjects={setSearchProjects}
            onChangeFlags={setSearchFlags}
            projectOptions={[...new Set(users.map((u) => u.projects).flat())]}
            flagOptions={[...new Set(users.map((u) => u.flags || []).flat())]}
          />
        </Grid>
        <Grid container className={classes.manage} item xs={3} alignItems='flex-end' justify='flex-start'>
          {selected.length > 0 ? (
            <Button onClick={openAddRemoveFlags}>Manage Flags</Button>
          ) : null}
          {addRemoveOpen && selected.length > 0 ? (
            <AddRemove selectedUsers={selected}/>
          ) : null}
        </Grid>
      </Grid>
      <Grid item xs={12}>
        <TableContainer component={Paper}>
          <Table className={classes.table} aria-label='user flags'>
            <TableHead>
              <TableRow>
                <TableCell padding='checkbox'>
                  <Checkbox
                    indeterminate={
                      selected.length > 0 &&
                      selected.length < filteredUsers.length
                    }
                    checked={
                      selected.length > 0 &&
                      selected.length === filteredUsers.length
                    }
                    onChange={onSelectAllClick}
                    inputProps={{'aria-label': 'select all users'}}
                  />
                </TableCell>
                <TableCell className={classes.header}>Name</TableCell>
                <TableCell className={classes.header}>Email</TableCell>
                <TableCell className={classes.header}>Projects</TableCell>
                <TableCell className={classes.header}>Flags</TableCell>
              </TableRow>
            </TableHead>
            <TableBody>
              {filteredUsers.map((user) => (
                <UserRow
                  user={user}
                  key={user.email}
                  onClick={onClickUser}
                  isSelected={isSelected(user)}
                />
              ))}
            </TableBody>
          </Table>
        </TableContainer>
      </Grid>
    </Grid>
  );
};

export default UserTable;
