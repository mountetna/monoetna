import React, {useState} from 'react';

import Menu from '@material-ui/core/Menu';
import MenuItem from '@material-ui/core/MenuItem';
import IconButton from '@material-ui/core/IconButton';
import MoreVertIcon from '@material-ui/icons/MoreVert';

export default function FigureMenu({
  figureId,
  onCopy,
  onRename,
  onRemove
}: {
  figureId: number;
  onCopy: () => void;
  onRename: () => void;
  onRemove: () => void;
}) {
  const [menuAnchor, setMenuAnchor] = useState(
    null as HTMLButtonElement | null
  );

  const handleClose = () => {
    setMenuAnchor(null);
  };

  return (
    <>
      <IconButton
        onClick={(e) => {
          e.stopPropagation();
          setMenuAnchor(e.currentTarget);
        }}
      >
        <MoreVertIcon />
      </IconButton>
      <Menu
        id={`figure-menu-${figureId}`}
        open={Boolean(menuAnchor)}
        anchorEl={menuAnchor}
        onClose={handleClose}
      >
        <MenuItem onClick={onCopy}>Copy</MenuItem>
        <MenuItem onClick={onRename}>Rename</MenuItem>
        <MenuItem onClick={onRemove}>Remove</MenuItem>
      </Menu>
    </>
  );
}
