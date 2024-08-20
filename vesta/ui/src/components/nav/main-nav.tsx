'use client'

import * as React from 'react'
import Container from '@mui/system/Container'
import Box from '@mui/system/Box'
import { Breakpoint, useTheme } from '@mui/material/styles';
import _ from 'lodash'
import { FocusTrap } from '@mui/base';
import { usePathname } from 'next/navigation';

import UCSFNav, { Classes as UCSFNavClasses } from './ucsf-nav'
import UCSFHomeLink from './ucsf-home-link';
import useIsStuck from '@/lib/utils/css';
import OverlayNav from './overlay-nav';
import NavBar, { Heights as NavBarHeights } from './nav-bar';
import { useBreakpoint } from '@/lib/utils/responsive';


const stickablePaths = ['/']


export function getMainNavHeight(breakpoint: Breakpoint): number {
    return NavBarHeights[breakpoint].condensed
}


export default function MainNav() {
    const theme = useTheme()

    const [overlayNavOpen, setOverlayNavOpen] = React.useState(false)
    const [dlNavFocus, setDlNavFocus] = React.useState(false)

    React.useEffect(() => {
        document.body.style.overflow = overlayNavOpen ? 'hidden' : 'visible'
    }, [overlayNavOpen])

    const breakpoint = useBreakpoint()
    const isDesktop = ['desktop', 'desktopLg'].includes(breakpoint)

    React.useEffect(() => {
        if (isDesktop && overlayNavOpen) {
            setOverlayNavOpen(false)
        }
    }, [isDesktop])

    const mainNavRef = React.createRef<HTMLElement>()
    const isStuck = useIsStuck(mainNavRef)
    const pathname = usePathname()
    const isStickablePath = stickablePaths.includes(pathname)
    const shouldAndIsStuck = !isStickablePath || isStuck

    const handleClickOverlayNavButton = () => {
        setOverlayNavOpen(!overlayNavOpen)
        setDlNavFocus(!overlayNavOpen)
    }

    // Compensate for anything above nav bar (i.e. UCSF nav)
    const [scrollDistanceToMainNav, setScrollDistanceToMainNav] = React.useState(0)

    React.useEffect(() => {
        const mainNavEl = mainNavRef.current
        if (!isDesktop && mainNavEl && !isStuck) {
            const setDistanceToMainNav = _.throttle(() => {
                setScrollDistanceToMainNav(
                    Math.max(
                        mainNavEl.getBoundingClientRect().top,
                        0,
                    )
                )
            }, 10)

            // Run once for initial page load
            setDistanceToMainNav()
            document.addEventListener('scroll', setDistanceToMainNav)

            return () => {
                document.removeEventListener('scroll', setDistanceToMainNav)
                setScrollDistanceToMainNav(0)
            }
        }
    }, [isDesktop, isStuck])

    return (
        <React.Fragment>
            {/* UCSF Nav */}
            <Box
                component='nav'
                aria-label='UCSF Main'
                sx={{
                    py: '10px',
                    bgcolor: 'utilityUCSFNavy.main',
                }}
            >
                <Container
                    sx={{
                        display: 'flex',
                        justifyContent: 'space-between',
                    }}
                >
                    <UCSFHomeLink />

                    <UCSFNav
                        sx={{
                            display: 'none',
                            [theme.breakpoints.up('desktop')]: {
                                display: 'flex',
                                ml: '29px',
                                gap: '29px',
                            },
                            [`.${UCSFNavClasses.link}`]: {
                                color: 'utilityHighlight.main',
                            },
                        }}
                    />
                </Container>
            </Box>

            <FocusTrap
                open={dlNavFocus}
            >
                {/* DL Nav */}
                <Box
                    ref={mainNavRef}
                    component='nav'
                    aria-label='Main'
                    zIndex='appBar'
                    sx={{
                        position: 'sticky',
                        top: 0,
                        bgcolor: 'utilityLowlight.main',
                        height: `${NavBarHeights[breakpoint].condensed}px`,
                        overflow: 'visible',
                        boxShadow: shouldAndIsStuck ? 'rgba(0, 0, 0, 0.15) 0px 0px 22px 0px' : 'rgba(0, 0, 0, 0) 0px 0px 22px 0px',
                        transition: theme.transitions.create(
                            'box-shadow',
                            {
                                duration: theme.transitions.duration.ease,
                                easing: theme.transitions.easing.ease,
                            },
                        ),
                        '& > *:first-child': {
                            [theme.breakpoints.up('desktop')]: {
                                display: 'none',
                            },
                        }
                    }}
                >
                    <OverlayNav
                        open={overlayNavOpen}
                        onClickNavLink={() => setOverlayNavOpen(false)}
                        sx={{
                            // Compensate for navbar height
                            pt: `${NavBarHeights[breakpoint][shouldAndIsStuck ? 'condensed' : 'default']}px`,
                            // Compensate for UCSF navbar visibile height
                            height: `calc(100vh - ${scrollDistanceToMainNav}px)`,
                            transition: theme.transitions.create(
                                'padding-top',
                                {
                                    duration: theme.transitions.duration.ease,
                                    easing: theme.transitions.easing.ease,
                                },
                            ),
                            '& > *:first-child': {
                                pt: '36px',
                            },
                        }}
                    />

                    <NavBar
                        variant={shouldAndIsStuck ? 'condensed' : 'default'}
                        sx={{
                            position: 'absolute',
                            top: 0,
                            width: '100%',
                        }}
                        onClickOverlayNavButton={handleClickOverlayNavButton}
                    />
                </Box>
            </FocusTrap>
        </React.Fragment >
    )
}