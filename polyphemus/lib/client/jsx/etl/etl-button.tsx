import React from 'react';
import Tooltip from '@material-ui/core/Typography';
import IconButton from '@material-ui/core/IconButton';

const EtlButton = ({children, mode, selected, onClick}:{children:React.ReactNode, mode:string, selected:string, onClick:Function}) => (
  <Tooltip title={mode}>
    <IconButton color={ mode == selected ? 'secondary' : 'default' } onClick={() => onClick(mode)} size='small' aria-label={mode}>
      {children}
    </IconButton>
  </Tooltip>
);

export default EtlButton;
