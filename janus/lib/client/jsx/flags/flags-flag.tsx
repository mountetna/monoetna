import React from 'react';

import {makeStyles} from '@material-ui/core/styles';
import Chip from '@material-ui/core/Chip';

const useStyles = makeStyles({});

const Flag = ({flag}: {flag: string}) => {
  const classes = useStyles();

  return <Chip label={flag} />;
};

export default Flag;
