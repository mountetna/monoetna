'use client'

import * as React from 'react';
import { SxProps, Typography, useTheme } from '@mui/material'


interface Props {
    sx?: SxProps
}


export default function TooltipContent({
    children,
    ...props
}: { children: React.ReactNode } & Props) {

    const {
        sx,
    } = props

    const theme = useTheme()

    return (
        <Typography
            component='div'
            variant='p2XSMono'
            sx={{
                m: '0.5em',
                px: '7px',
                py: '2px',
                borderRadius: '8px',
                bgcolor: 'blue.grade50',
                color: 'utilityWhite.main',
                ...(sx ? sx : {}),
            }}
        >
            {children}
        </Typography>
    )
}