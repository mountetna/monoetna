import * as React from 'react';
import ButtonBase from '@mui/material/ButtonBase';


export enum Classes {
    root = 'library-card-button'
}


export default function LibraryCardButton({
    isLoggedIn,
    textOverride,
    onClick,
}: {
    isLoggedIn: boolean,
    textOverride?: string,
    onClick: () => void,
}) {
    return (
        <ButtonBase
            className={Classes.root}
            onClick={onClick}
            sx={{
                width: 'fit-content',
                px: '16px',
                py: '8px',
                color: 'utilityHighlight.main',
                backgroundColor: 'blue.grade50',
                borderRadius: '10px',
                typography: 'pBodyMediumWt',
            }}
        >
            {textOverride !== undefined ? textOverride : isLoggedIn ? 'View your Library Card' : 'Get Access'}
        </ButtonBase>
    )
}