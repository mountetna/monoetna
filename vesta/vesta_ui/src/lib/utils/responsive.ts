'use client'

import * as React from 'react'
import _ from 'lodash'
import { Breakpoint, useMediaQuery, useTheme } from '@mui/material'


export function useWindowDimensions(triggerDelayMs: number = 100) {
    const [dimensions, setDimensions] = React.useState<number[] | null>(null)
    const [isResizing, setIsResizing] = React.useState<boolean>(false)

    React.useEffect(() => {
        const dimensionsHandler = _.throttle(() => {
            setDimensions([window.innerWidth, window.innerHeight])
        }, triggerDelayMs)

        const resizingTrueHandler = _.throttle(() => {
            setIsResizing(true)
        }, triggerDelayMs)

        const resizingFalseHandler = _.debounce(() => {
            resizingTrueHandler.cancel()
            setIsResizing(false)
        }, triggerDelayMs)

        const handlers = [
            dimensionsHandler,
            resizingTrueHandler,
            resizingFalseHandler,
        ]

        for (const handler of handlers) {
            window.addEventListener('resize', handler)
            // trigger handlers on react mount
            handler()
        }

        return () => {
            for (const handler of handlers) {
                window.removeEventListener('resize', handler)
            }
        }
    }, [triggerDelayMs])

    return { dimensions, isResizing }
}

export function useBreakpoint(): Breakpoint {
    const theme = useTheme()
    const isMobile = useMediaQuery(theme.breakpoints.down('tablet'))
    const isTablet = useMediaQuery(theme.breakpoints.between('tablet', 'desktop'))
    const isDesktop = useMediaQuery(theme.breakpoints.between('desktop', 'desktopLg'))
    const isDesktopLg = useMediaQuery(theme.breakpoints.up('desktopLg'))

    function getBreakpoint(): Breakpoint {
        if (isMobile) return 'mobile'
        else if (isTablet) return 'tablet'
        else if (isDesktop) return 'desktop'
        else return 'desktopLg'
    }

    const [breakpoint, setBreakpoint] = React.useState(getBreakpoint())

    React.useEffect(() => {
        setBreakpoint(getBreakpoint())
    }, [isMobile, isTablet, isDesktop, isDesktopLg])

    return breakpoint
}