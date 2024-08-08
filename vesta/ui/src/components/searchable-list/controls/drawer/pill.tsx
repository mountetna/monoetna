'use client'

import * as React from 'react'
import Typography from '@mui/material/Typography';
import ButtonBase from '@mui/material/ButtonBase';
import { useTheme } from '@mui/material';


export enum Classes {
    base = 'drawer-pill',
    active = 'drawer-pill-active',
}


export default function DrawerPill({
    label,
    active = false,
    onClick,
    variant,
}: {
    label: string,
    active?: boolean,
    onClick?: () => void,
    variant: 'yellow' | 'teal',
}) {
    const theme = useTheme()

    const classes = [Classes.base]
    if (active) {
        classes.push(Classes.active)
    }

    return (
        <ButtonBase
            className={classes.join(' ')}
            onClick={() => onClick && onClick()}
            sx={{
                px: '12px',
                py: '2px',
                color: 'ground.grade10',
                bgcolor: 'utilityHighlight.main',
                borderRadius: '40px',
                border: '1px solid transparent',
                [`&.${Classes.active}`]: {
                    bgcolor: variant === 'yellow' ? 'yellow.grade50' : 'teal.grade100',
                },
                transition: theme.transitions.create(
                    'all',
                    {
                        easing: theme.transitions.easing.ease,
                        duration: theme.transitions.duration.ease,
                    },
                ),
            }}
        >
            <Typography variant='pMedium'>
                {label}
            </Typography>
        </ButtonBase>
    )
}