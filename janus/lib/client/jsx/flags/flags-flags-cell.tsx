import React, {useState, useEffect, useCallback} from 'react';
import * as _ from 'lodash';

import {makeStyles} from '@material-ui/core/styles';
import AddIcon from '@material-ui/icons/Add';
import IconButton from '@material-ui/core/IconButton';
import Tooltip from '@material-ui/core/Tooltip';

import Flag from './flags-flag';

const useStyles = makeStyles((theme) => ({
  margin: {
    margin: theme.spacing(3)
  }
}));

const FlagsCell = ({
  flags,
  onChange
}: {
  flags: string[] | null;
  onChange: (flags: string[]) => void;
}) => {
  return (
    <React.Fragment>
      {flags?.map((flag) => (
        <Flag flag={flag} />
      ))}
      {/* <Tooltip title='Add'>
        <IconButton aria-label='add flag'>
          <AddIcon />
        </IconButton>
      </Tooltip> */}
    </React.Fragment>
  );
};

export default FlagsCell;
