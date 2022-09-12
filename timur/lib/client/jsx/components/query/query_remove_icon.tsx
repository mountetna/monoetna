import React from 'react';

import IconButton from '@material-ui/core/IconButton';
import ClearIcon from '@material-ui/icons/Clear';
import Tooltip from '@material-ui/core/Tooltip';

const RemoveIcon = React.memo(
  ({
    showRemoveIcon,
    onClick,
    label
  }: {
    showRemoveIcon: boolean;
    onClick: () => void;
    label: string;
  }) => {
    if (!showRemoveIcon) return null;

    return (
      <Tooltip title={`Remove ${label}`} aria-label={`remove ${label}`}>
        <IconButton aria-label={`remove ${label}`} onClick={onClick}>
          <ClearIcon color='action' />
        </IconButton>
      </Tooltip>
    );
  }
);

export default RemoveIcon;
