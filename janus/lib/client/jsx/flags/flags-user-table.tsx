import React, {useState, useEffect, useCallback} from 'react';
import Grid from '@material-ui/core/Grid';
import {makeStyles} from '@material-ui/core/styles';
import TextField from '@material-ui/core/TextField';
import InputAdornment from '@material-ui/core/InputAdornment';
import Search from '@material-ui/icons/Search';

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
  const classes = useStyles();

  return (
    <TableContainer component={Paper}>
      <Table className={classes.table} aria-label='user flags'>
        <TableHead>
          <TableRow>
            <TableCell className={classes.header}>Name</TableCell>
            <TableCell className={classes.header}>Email</TableCell>
            <TableCell className={classes.header}>Projects</TableCell>
            <TableCell className={classes.header}>Flags</TableCell>
            <TableCell className={classes.header}></TableCell>
          </TableRow>
        </TableHead>
        <TableBody>
          {users.map((user) => (
            <UserRow user={user} key={user.email} />
          ))}
        </TableBody>
      </Table>
    </TableContainer>
  );
};

export default UserTable;
