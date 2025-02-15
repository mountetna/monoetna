'use client'

import * as React from 'react';
import { Box, SxProps, Typography, BoxProps } from '@mui/material'


interface Props extends BoxProps {
    variant: 'simple' | 'expanded'
}


export default function TooltipContent({
    children,
    ...props
}: { children: React.ReactNode } & Props) {

    const {
        sx,
        variant,
    } = props

    switch (variant) {
        case 'simple':

            return (
                <Box
                    {...props}
                    sx={{
                        p: '0.5em',
                        ...(sx ? sx : {}),
                    }}
                >
                    <Box
                        sx={{
                            px: '7px',
                            py: '2px',
                            borderRadius: '8px',
                            bgcolor: 'blue.grade50',
                            color: 'utilityWhite.main',
                        }}
                    >
                        <Typography
                            variant='p2XSMono'
                            sx={{
                                color: 'inherit',
                            }}
                        >
                            {children}
                        </Typography>
                    </Box>
                </Box>
            )

        case 'expanded':

            return (
                <Box
                    {...props}
                    sx={{
                        p: '10px',
                        ...(sx ? sx : {}),
                    }}
                >
                    <Box
                        sx={{
                            display: 'flex',
                            flexDirection: 'column',
                            p: '10px',
                            gap: '10px',
                            borderRadius: '16px',
                            border: '1px solid',
                            borderColor: 'utilityHighlight.main',
                            bgcolor: 'utilityWhite.main',
                        }}
                    >
                        {children}
                    </Box>
                </Box>
            )
    }
}