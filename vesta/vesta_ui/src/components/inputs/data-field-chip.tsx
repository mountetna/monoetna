'use client'

import * as React from 'react'
import Typography from '@mui/material/Typography';
import ButtonBase from '@mui/material/ButtonBase';
import { SxProps, useTheme } from '@mui/material';
import { TypographyVariant } from '@/lib/utils/types';


export default function Button({
    label,
    onClick,
    sizeVariant,
    strokeVariant,
    typographyVariant = 'pSubtitle',
    sx,
    disabled = false,
}: {
    label: string,
    onClick: () => void,
    sizeVariant: 'small' | 'large',
    strokeVariant: 'filled' | 'stroked',
    typographyVariant: TypographyVariant,
    sx?: SxProps,
    disabled?: boolean,
}) {
    const theme = useTheme()

    return (
        <ButtonBase
            onClick={() => onClick()}
            sx={{
                width: '328.57px',
                height: '65.71px',
                px: '16.4px',
                py: '8.21px',
                borderRadius: '8.21px',
                border: `2px solid ${strokeVariant === 'filled' ? 'transparent' : theme.palette.ground.grade10}`,
                bgcolor: strokeVariant === 'filled' ? 'ground.grade10' : 'utilityWhite.main',
                color: strokeVariant === 'filled' ? '#F5F5F5' : 'ground.grade10',
                '&:hover, &:focus': {
                    bgcolor: strokeVariant === 'filled' ? 'ground.grade25' : 'ground.grade10',
                    color: strokeVariant === 'filled' ? 'ground.grade10' : '#F5F5F5',
                },
                transition: theme.transitions.create(
                    'all',
                    {
                        easing: theme.transitions.easing.ease,
                        duration: theme.transitions.duration.ease,
                    },
                ),
                ...(sx ? sx : {}),
            }}
            disabled={disabled}
        >
            <Typography variant={typographyVariant}>
                {label}
            </Typography>
        </ButtonBase>
    )
}
