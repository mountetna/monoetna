'use client'

import * as React from 'react';
import MUILink, { LinkProps } from '@mui/material/Link';
import NextLink from 'next/link'
import Tooltip from '../tooltip/tooltip';
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
        id,
    } = props

    // Manage getting full href
    // (using `id` since MUILink doesn't support `ref`)
    let elID = React.useId()
    elID = id !== undefined ? id : elID

    const [fullHref, setFullHref] = React.useState(href)

    React.useEffect(() => {
        const anchorEl = document.getElementById(elID) as HTMLAnchorElement | null
        if (!anchorEl) return

        setFullHref(anchorEl.href)
    }, [])

    const link = (
        <MUILink
            {...props}
            id={elID}
            component={NextLink}
            underline='none'
        >
            {children}
        </MUILink>
    )

    return (
        tooltip ? (
            <Tooltip
                title={tooltipContent ? tooltipContent : (
                    <TooltipContent
                        variant='simple'
                    >
                        {fullHref}
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