import React, {useState, useEffect, useCallback} from 'react';
import * as _ from 'lodash';

import {makeStyles} from '@material-ui/core/styles';
import {UserFlagsInterface} from '../models/user_models';
import {updateUserFlags} from '../api/janus_api';
import {TableCell} from '@material-ui/core';
import TableRow from '@material-ui/core/TableRow';

import FlagsCell from './flags-flags-cell';

const useStyles = makeStyles((theme) => ({
  cell: {
    maxWidth: 300,
    overflow: 'hidden',
    minWidth: 200,
    textOverflow: 'ellipsis'
  }
}));

const UserRow = ({user}: {user: UserFlagsInterface}) => {
  const [updatedFlags, setUpdatedFlags] = useState([] as string[] | null);
  const [allowSave, setAllowSave] = useState(false as boolean);

  const classes = useStyles();

  useEffect(() => {
    setUpdatedFlags(user.flags);
  }, [user.flags]);

  useEffect(() => {
    setAllowSave(!_.isEqual(updatedFlags, user.flags));
  }, [updatedFlags]);

  return (
    <TableRow>
      <TableCell className={classes.cell}>{user.name}</TableCell>
      <TableCell className={classes.cell}>{user.email}</TableCell>
      <TableCell className={classes.cell}>{user.projects.join(',')}</TableCell>
      <TableCell className={classes.cell}>
        <FlagsCell flags={updatedFlags} onChange={setUpdatedFlags} />
      </TableCell>
      <TableCell className={classes.cell}>
        {allowSave ? <div> Save buttn </div> : null}
      </TableCell>
    </TableRow>
  );
};

export default UserRow;
