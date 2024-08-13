'use client'

import * as React from 'react';
import MUILink, { LinkProps } from '@mui/material/Link';
import Box from '@mui/system/Box'
import NextLink from 'next/link'
import { Tooltip as MUITooltip } from '@mui/material';


interface Props extends Omit<LinkProps, 'component' | 'children' | 'underline' | 'href'> {
    href: string
    tooltip?: boolean
}


export default function Link({
    children,
    ...props
}: { children: React.ReactNode } & Props) {
    const link = (
        <MUILink
            component={NextLink}
            underline='none'
            {...props}
        >
            {children}
        </MUILink>
    )

    return (
        props.tooltip || props.tooltip === undefined ? (
            <MUITooltip
                title={props.href}
                followCursor
            >
                {link}
            </MUITooltip>
        ) : (
            link
        )
    )
}