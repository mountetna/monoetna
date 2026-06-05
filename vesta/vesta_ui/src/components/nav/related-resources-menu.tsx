import * as React from 'react'
import Menu from '@mui/material/Menu';
import MenuItem from '@mui/material/MenuItem';

import { TypographyVariant } from '@/lib/utils/types'
import NavLink from './nav-link';

export default function RelatedResourcesMenu({
    typography = 'pBody',
    anchorEl,
    onClick,
    onClose
}: {
    typography?: TypographyVariant,
    anchorEl: any,
    onClick: (event: React.MouseEvent<HTMLAnchorElement>) => void,
    onClose: () => void,
}) {
    const open = Boolean(anchorEl);

    return <Menu anchorEl={anchorEl} open={open} onClose={onClose} >
      <MenuItem onClick={onClose} sx={{ '&:hover': { backgroundColor: 'white' } }}>
        <NavLink
          text='Legacy Version'
          href={'https://datalibraryarchive.ucsf.edu/'}
          onClick={onClick}
          typography={typography}
        />
      </MenuItem>
      <MenuItem onClick={onClose} sx={{ '&:hover': { backgroundColor: 'white' } }}>
        <NavLink
          text='Biobanks'
          href={'https://ocr.ucsf.edu/biobanks'}
          onClick={onClick}
          typography={typography}
        />
      </MenuItem>
      <MenuItem onClick={onClose} sx={{ '&:hover': { backgroundColor: 'white' } }}>
        <NavLink
          text='ImmunoX Data Commons'
          href={'https://immunox.notion.site/ImmunoX-Data-Commons-118f8b14a9b780478de8c1e084d577bc?pvs=74'}
          onClick={onClick}
          typography={typography}
        />
      </MenuItem>
      <MenuItem onClick={onClose} sx={{ '&:hover': { backgroundColor: 'white' } }}>
        <NavLink
          text='Analysis Tools'
          href={'https://ocr.ucsf.edu/analysis-tools'}
          onClick={onClick}
          typography={typography}
        />
      </MenuItem>
    </Menu>
}
