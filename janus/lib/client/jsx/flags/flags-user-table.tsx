import React, {useState, useEffect, useCallback} from 'react';

import {makeStyles} from '@material-ui/core/styles';
import Checkbox from '@material-ui/core/Checkbox';
import Grid from '@material-ui/core/Grid';
import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableContainer from '@material-ui/core/TableContainer';
import TableHead from '@material-ui/core/TableHead';
import TableRow from '@material-ui/core/TableRow';
import Paper from '@material-ui/core/Paper';

import {UserFlagsInterface} from '../types/janus_types';
import UserRow from './flags-user-row';
import TableControls from './flags-table-controls';
import AddRemove from './flags-add-remove';

const useStyles = makeStyles((theme) => ({
  table: {
    minWidth: 650
  },
  header: {
    fontWeight: 'bolder'
  },
  margin: {
    margin: theme.spacing(3)
  }
}));

const UserTable = ({
  users,
  onChange
}: {
  users: UserFlagsInterface[];
  onChange: () => void;
}) => {
  const [filteredUsers, setFilteredUsers] = useState(
    [] as UserFlagsInterface[]
  );
  const [searchTerm, setSearchTerm] = useState('');
  const [searchProjects, setSearchProjects] = useState([] as string[]);
  const [searchFlags, setSearchFlags] = useState([] as string[]);
  const [selected, setSelected] = useState([] as UserFlagsInterface[]);
  const classes = useStyles();

  useEffect(() => {
    setFilteredUsers(users);
  }, [users]);

  useEffect(() => {
    // (searchTerm across user.name || user.email) &&
    //  (searchProjects OR'd) && (searchFlags OR'd)
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
  }, [searchTerm, searchProjects, searchFlags]);

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

  return (
    <Grid container xs={12} direction='column' className={classes.margin}>
      <Grid container item>
        <Grid item xs={9}>
          <TableControls
            onChangeSearch={setSearchTerm}
            onChangeProjects={setSearchProjects}
            onChangeFlags={setSearchFlags}
            projectOptions={[...new Set(users.map((u) => u.projects).flat())]}
            flagOptions={[...new Set(users.map((u) => u.flags || []).flat())]}
          />
        </Grid>
        <Grid item xs={3}>
          {selected.length > 0 ? (
            <AddRemove
              selectedUsers={selected}
              onUpdateComplete={() => {
                setSelected([]);
                onChange();
              }}
            />
          ) : null}
        </Grid>
      </Grid>
      <Grid item>
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
                <TableCell className={classes.header}></TableCell>
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
