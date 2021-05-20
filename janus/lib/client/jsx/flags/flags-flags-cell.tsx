import React, {useState, useEffect, useCallback} from 'react';

import {makeStyles} from '@material-ui/core/styles';

import Flag from './flags-flag';

const useStyles = makeStyles((theme) => ({
  margin: {
    margin: theme.spacing(3)
  }
}));

const FlagsCell = ({flags}: {flags: string[] | null}) => {
  return (
    <React.Fragment>
      {flags?.map((flag) => (
        <Flag flag={flag} />
      ))}
    </React.Fragment>
  );
};

export default FlagsCell;
