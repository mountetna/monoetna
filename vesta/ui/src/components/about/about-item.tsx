import * as React from 'react'
import Box from '@mui/system/Box'
import Container from '@mui/system/Container'
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
            sx={{
                ...containerPadding,
                width: '100vw',
                maxWidth: 'desktopLg',
                backgroundColor: 'black',
            }}
        >
            {header}
        </Box>
    )
}