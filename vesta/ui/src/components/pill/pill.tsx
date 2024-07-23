'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { useTheme } from '@mui/material';
import Image, { StaticImageData } from 'next/image';

import { TypographyVariant } from '@/lib/utils/types';


export default function Pill({
    label,
    typographyVariant = 'pMedium',
    icon,
    iconAlt,
    iconPosition = 'before',
    variant,
}: {
    label: string,
    typographyVariant: TypographyVariant,
    icon?: StaticImageData,
    iconAlt?: string,
    iconPosition?: 'before' | 'after',
    variant: 'filled' | 'stroked',
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

    return (
        <Box
            className={`pill pill-${variant}`}
            sx={{
                display: 'flex',
                flexDirection: 'row',
                gap: '6px',
                px: '10px',
                py: '2px',
                borderRadius: '20px',
                alignItems: 'center',
                textAlign: 'center',
                '&.pill-filled': {
                    background: theme.palette.utilityWhite.main,
                },
                '&.pill-stroked': {
                    border: `1px solid ${theme.palette.ground.grade10}`,
                    background: 'transparent',
                },
            }}
        >
            {iconPosition === 'before' && iconEl}

            <Typography
                variant={typographyVariant}
            >
                {label}
            </Typography>

            {iconPosition === 'after' && iconEl}
        </Box>
    )
}