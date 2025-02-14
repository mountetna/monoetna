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
    counter = 'pill-counter',
}


export default function Pill({
    label,
    typographyVariant = 'pMedium',
    icon,
    iconAlt,
    iconPosition = 'before',
    // TODO: make CounterPill its own class
    count,
    variant,
    sx = {},
}: {
    label: string,
    typographyVariant: TypographyVariant,
    icon?: StaticImageData,
    iconAlt?: string,
    iconPosition?: 'before' | 'after',
    count?: number,
    variant: 'filled' | 'stroked' | 'counter',
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

    const classes = [Classes.base, Classes[variant]]

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
                [`&.${Classes.counter}`]: {
                    border: `1px solid ${theme.palette.ground.grade10}`,
                    background: 'transparent',
                    py: '3px',
                    pl: '3px',
                },
                ...sx,
            }}
        >
            {variant === 'counter' && (
                <Typography
                    variant={typographyVariant}
                    sx={{
                        display: 'inline-flex',
                        color: 'utilityWhite.main',
                        bgcolor: 'utilityLowlight.main',
                        borderRadius: '40px',
                        px: '7px',
                    }}
                >
                    {count}
                </Typography>
            )}

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