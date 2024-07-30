import * as React from 'react';
import Box from '@mui/system/Box'
import ButtonBase from '@mui/material/ButtonBase';


export enum Classes {
    root = 'library-card-button'
}


export default function LibraryCardButton() {
    return (
        <ButtonBase
            className={Classes.root}
            sx={{
                px: '16px',
                py: '8px',
                color: 'utilityHighlight.main',
                backgroundColor: 'blue.grade50',
                borderRadius: '10px',
                typography: 'pBodyMediumWt',
            }}
        >
            LIB CARD / GET ACCESS
        </ButtonBase>
    )
}