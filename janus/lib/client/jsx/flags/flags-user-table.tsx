import React, {useState, useEffect, useCallback} from 'react';

import {makeStyles} from '@material-ui/core/styles';
import Checkbox from '@material-ui/core/Checkbox';
import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableContainer from '@material-ui/core/TableContainer';
import TableHead from '@material-ui/core/TableHead';
import TableRow from '@material-ui/core/TableRow';
import Paper from '@material-ui/core/Paper';

import {UserFlagsInterface} from '../models/user_models';
import UserRow from './flags-user-row';

const useStyles = makeStyles({
  table: {
    minWidth: 650
  },
  header: {
    fontWeight: 'bolder'
  }
});

const UserTable = ({users}: {users: UserFlagsInterface[]}) => {
  const [selected, setSelected] = useState([] as UserFlagsInterface[]);
  const classes = useStyles();

  function onSelectAllClick() {
    setSelected(selected.length > 0 ? [] : users);
  }

  function onClickUser(event: any, user: UserFlagsInterface) {
    if (event.target.checked) {
      setSelected([...selected].concat([user]));
    } else {
      setSelected([...selected].filter((u) => u.email !== user.email));
    }
  }

  function isSelected(user: UserFlagsInterface) {
    return selected.filter((u) => u.email === user.email).length > 0;
  }

  return (
    <TableContainer component={Paper}>
      <Table className={classes.table} aria-label='user flags'>
        <TableHead>
          <TableRow>
            <TableCell padding='checkbox'>
              <Checkbox
                indeterminate={
                  selected.length > 0 && selected.length < users.length
                }
                checked={
                  selected.length > 0 && selected.length === users.length
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
          {users.map((user) => (
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
  );
};

export default UserTable;
