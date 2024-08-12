'use client'

import * as React from 'react'
import Typography from '@mui/material/Typography';
import ButtonBase from '@mui/material/ButtonBase';
import { useTheme } from '@mui/material';
import Pill, { Classes as PillClasses } from '@/components/pill/pill';


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
                [`&.${Classes.active} .${PillClasses.filled}`]: {
                    bgcolor: variant === 'yellow' ? 'yellow.grade50' : 'teal.grade100',
                },
            }}
        >
            <Pill
                label={label}
                typographyVariant='pMedium'
                variant='filled'
                sx={{
                    color: 'ground.grade10',
                    borderRadius: '40px',
                    border: '1px solid transparent',
                    px: '12px',
                    py: '2px',
                    [`&.${PillClasses.filled}`]: {
                        bgcolor: 'utilityHighlight.main',
                    },
                    transition: theme.transitions.create(
                        'all',
                        {
                            easing: theme.transitions.easing.ease,
                            duration: theme.transitions.duration.ease,
                        },
                    ),
                }}
            />
        </ButtonBase>
    )
}