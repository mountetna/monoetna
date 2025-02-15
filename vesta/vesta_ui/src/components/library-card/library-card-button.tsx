import * as React from 'react';
import ButtonBase from '@mui/material/ButtonBase';
import { SxProps, Typography } from '@mui/material';

import Link from '../link/link';


export enum Classes {
    root = 'library-card-button'
}


export default function LibraryCardButton({
    isLoggedIn,
    accessUrl,
    textOverride,
    onClick,
    disabled = false,
}: {
    isLoggedIn: boolean,
    accessUrl: string,
    textOverride?: string,
    onClick: () => void,
    disabled?: boolean,
}) {
    const textEl = (
        <Typography variant='pBodyMediumWt'>
            {textOverride !== undefined ? textOverride : isLoggedIn ? 'View your Library Card' : 'Get Access'}
        </Typography>
    )

    const sx: SxProps = {
        width: 'fit-content',
        px: '16px',
        py: '8px',
        color: 'utilityHighlight.main',
        backgroundColor: 'blue.grade50',
        borderRadius: '10px',
    }

    return isLoggedIn ? (
        <ButtonBase
            tabIndex={0}
            className={Classes.root}
            onClick={onClick}
            sx={sx}
            disabled={disabled}
        >
            {textEl}
        </ButtonBase>
    ) : (
        <Link
            href={accessUrl}
            sx={sx}
        >
            {textEl}
        </Link>
    )
}