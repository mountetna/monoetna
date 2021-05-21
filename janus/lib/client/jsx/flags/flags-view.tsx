import React from 'react';
import ScopedCssBaseline from '@material-ui/core/ScopedCssBaseline';
import {makeStyles} from '@material-ui/core/styles';

import {FlagsProvider} from './flags-context';

import UserTable from './flags-user-table';

const useStyles = makeStyles((theme) => ({
  margin: {
    margin: theme.spacing(3)
  }
}));

const FlagsView = () => {
  const classes = useStyles();

  return (
    <ScopedCssBaseline>
      <FlagsProvider>
        <UserTable />
      </FlagsProvider>
    </ScopedCssBaseline>
  );
};

export default FlagsView;
