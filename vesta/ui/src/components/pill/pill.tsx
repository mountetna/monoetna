'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { SxProps, useTheme } from '@mui/material';
import Image, { StaticImageData } from 'next/image';

import { TypographyVariant } from '@/lib/utils/types';


export enum Classes {
    base = 'pill',
    filled = 'pill-filled',
    stroked = 'pill-stroked',
}


export default function Pill({
    label,
    typographyVariant = 'pMedium',
    icon,
    iconAlt,
    iconPosition = 'before',
    variant,
    sx = {},
}: {
    label: string,
    typographyVariant: TypographyVariant,
    icon?: StaticImageData,
    iconAlt?: string,
    iconPosition?: 'before' | 'after',
    variant: 'filled' | 'stroked',
    sx?: SxProps,
}) {
    const theme = useTheme()

    let iconEl
    if (icon) {
        iconEl = (
            <Image
                src={icon}
                alt={iconAlt ? iconAlt : `Image for ${label}`}
            />
        )
    }

    const classes = [Classes.base]
    switch (variant) {
        case 'filled':
            classes.push(Classes.filled)
            break
        case 'stroked':
            classes.push(Classes.stroked)
            break
    }

    return (
        <Box
            className={classes.join(' ')}
            sx={{
                position: 'relative',
                display: 'flex',
                flexDirection: 'row',
                gap: '6px',
                px: '10px',
                py: '2px',
                borderRadius: '20px',
                alignItems: 'center',
                textAlign: 'center',
                [`&.${Classes.filled}`]: {
                    background: theme.palette.utilityWhite.main,
                },
                [`&.${Classes.stroked}`]: {
                    border: `1px solid ${theme.palette.ground.grade10}`,
                    background: 'transparent',
                },
                ...sx,
            }}
        >
            {iconPosition === 'before' && iconEl}

            <Typography
                variant={typographyVariant}
                sx={{
                    display: 'inline-flex',
                    color: 'ground.grade10',
                }}
            >
                {label}
            </Typography>

            {iconPosition === 'after' && iconEl}
        </Box>
    )
}