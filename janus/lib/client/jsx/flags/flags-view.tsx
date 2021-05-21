import React from 'react';
import ScopedCssBaseline from '@material-ui/core/ScopedCssBaseline';

import {FlagsProvider} from './flags-context';

import UserTable from './flags-user-table';

const FlagsView = () => {
  return (
    <ScopedCssBaseline>
      <FlagsProvider>
        <UserTable />
      </FlagsProvider>
    </ScopedCssBaseline>
  );
};

export default FlagsView;
