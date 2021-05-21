import React from 'react';

import Chip from '@material-ui/core/Chip';

const FlagsCell = ({flags}: {flags: string[] | null}) => {
  return (
    <React.Fragment>
      {flags?.map((flag, index) => (
        <Chip label={flag} key={index} />
      ))}
    </React.Fragment>
  );
};

export default FlagsCell;
