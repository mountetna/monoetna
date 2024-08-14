'use client'

import * as React from 'react';
import MUILink, { LinkProps } from '@mui/material/Link';
import NextLink from 'next/link'
import Tooltip from '../tooltip/tooltip';
import { useTheme } from '@mui/material';
import TooltipContent from '../tooltip/tooltip-content';


interface Props extends Omit<LinkProps, 'component' | 'children' | 'underline' | 'href'> {
    href: string
    tooltip?: boolean
    tooltipContent?: React.ReactNode
}


export default function Link({
    children,
    ...props
}: { children: React.ReactNode } & Props) {

    const {
        href,
        tooltip,
        tooltipContent,
    } = props

    const theme = useTheme()

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
        tooltip || tooltip === undefined ? (
            <Tooltip
                title={tooltipContent ? tooltipContent : (
                    <TooltipContent>
                        {href}
                    </TooltipContent>
                )}
                followCursor
                placement='bottom-start'
            >
                {link}
            </Tooltip>
        ) : (
            link
        )
    )
}