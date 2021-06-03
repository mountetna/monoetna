import React from 'react';

import {FlagsProvider} from './flags-context';

import UserTable from './flags-user-table';

const FlagsView = () => {
  return (
    <FlagsProvider>
      <UserTable />
    </FlagsProvider>
  );
};

export default FlagsView;
