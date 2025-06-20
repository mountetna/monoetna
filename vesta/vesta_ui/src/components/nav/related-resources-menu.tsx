import * as React from 'react'
import Menu from '@mui/material/Menu';
import MenuItem from '@mui/material/MenuItem';

import { TypographyVariant } from '@/lib/utils/types'
import NavLink from './nav-link';

export default function RelatedResourcesMenu({
    linkTypography = 'pBody',
    anchorEl,
    onClick,
    onClose
}: {
    linkTypography?: TypographyVariant,
    anchorEl: any,
    onClick: () => void,
    onClose: () => void,
}) {
    const open = Boolean(anchorEl);

    return <Menu anchorEl={anchorEl} open={open} onClose={onClose} >
      <MenuItem onClick={onClose} sx={{ '&:hover': { backgroundColor: 'white' } }}>
        <NavLink
          text='Legacy Version'
          href={'https://datalibraryarchive.ucsf.edu/'}
          onClick={onClick}
          typography={linkTypography}
        />
      </MenuItem>
      <MenuItem onClick={onClose} sx={{ '&:hover': { backgroundColor: 'white' } }}>
        <NavLink
          text='QuIPI'
          href={'https://quipi.org/'}
          onClick={onClick}
          typography={linkTypography}
        />
      </MenuItem>
    </Menu>
}
