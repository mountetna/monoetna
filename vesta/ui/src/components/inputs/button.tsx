'use client'

import * as React from 'react'
import Typography from '@mui/material/Typography';
import ButtonBase from '@mui/material/ButtonBase';
import { SxProps, useTheme } from '@mui/material';
import { TypographyVariant } from '@/lib/utils/types';


export default function Button({
    label,
    onClick,
    variant,
    typographyVariant = 'pSubtitle',
    sx,
    disabled = false,
}: {
    label: string,
    onClick: () => void,
    variant: 'small' | 'large',
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
                bgcolor: 'ground.grade10',
                color: '#F5F5F5',
                '&:hover, &:focus': {
                    bgcolor: 'ground.grade25',
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