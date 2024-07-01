import * as React from 'react';
import Box from '@mui/system/Box'
import ButtonBase from '@mui/material/ButtonBase';


// sx={{
//     display: 'inline-block',
//     px: '16px',
//     py: '8px',
//     borderRadius: '40px',
//     color: 'utilityUCSFNavy.main',
//     backgroundColor: 'utilityUCSFLightNavy.main',
// }}


export default function LibraryCardButton() {
    return (
        <ButtonBase
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