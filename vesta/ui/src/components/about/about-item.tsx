'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import { containerPadding } from '@/theme'


export interface Link {
    header?: string
    blurb?: string
    href: string
    label: string
}

export default function AboutItem({
    header,
    body,
    link,
    imageSrc,
}: {
    header: string,
    body: string,
    link?: Link,
    imageSrc: string,
}) {
    return (
        <Box
            className='about-item'
            sx={{
                width: '100vw',
                maxWidth: 'desktopLg',
            }}
        >
            <Box
                sx={{
                    ...containerPadding,
                }}
            >
                <Box
                    sx={(theme) => ({
                        bgcolor: 'utilityWhite.main',
                        borderRadius: '30px',
                        px: '8px',
                        py: '8px',
                        [theme.breakpoints.up('tablet')]: {
                            px: '30px',
                            py: '30px',
                        },
                        [theme.breakpoints.up('desktop')]: {
                            px: '30px',
                            py: '30px',
                        },
                    })}
                >
                    {header}
                </Box>
            </Box>
        </Box>
    )
}