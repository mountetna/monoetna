import * as React from 'react';
import Box from '@mui/system/Box'

import LibraryCardButton from './library-card-button';


export enum Classes {
    root = 'library-card-tray'
}


export default function LibraryCardTray() {
    return (
        <Box
            className={Classes.root}
        >
            <LibraryCardButton />
        </Box>
    )
}