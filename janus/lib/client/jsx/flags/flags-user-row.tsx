import React, {useState, useEffect, useCallback} from 'react';
import * as _ from 'lodash';
import Checkbox from '@material-ui/core/Checkbox';

import {makeStyles} from '@material-ui/core/styles';
import {UserFlagsInterface} from '../types/janus_types';
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

const UserRow = ({
  user,
  isSelected,
  onClick
}: {
  user: UserFlagsInterface;
  isSelected: boolean;
  onClick: (user: UserFlagsInterface) => void;
}) => {
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
    <TableRow
      hover
      onClick={() => onClick(user)}
      role='checkbox'
      aria-checked={isSelected}
      tabIndex={-1}
      selected={isSelected}
    >
      <TableCell padding='checkbox'>
        <Checkbox
          checked={isSelected}
          inputProps={{'aria-labelledby': user.email}}
        />
      </TableCell>
      <TableCell className={classes.cell}>{user.name}</TableCell>
      <TableCell className={classes.cell}>{user.email}</TableCell>
      <TableCell className={classes.cell}>{user.projects.join(',')}</TableCell>
      <TableCell className={classes.cell}>
        <FlagsCell flags={updatedFlags} />
      </TableCell>
    </TableRow>
  );
};

export default UserRow;
