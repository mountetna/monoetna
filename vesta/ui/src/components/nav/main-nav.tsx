'use client'

import * as React from 'react'
import Container from '@mui/system/Container'
import Box from '@mui/system/Box'
import { useMediaQuery } from '@mui/material';
import { useTheme } from '@mui/material/styles';
import _ from 'lodash'

import UCSFNav from './ucsf-nav'
import UCSFHomeLink from './ucsf-home-link';
import useIsStuck from '@/lib/utils/css';
import MobileNav from './mobile-nav';
import NavBar, { Heights as NavBarHeights, Classes as NavBarClasses } from './nav-bar';


export default function MainNav() {
    const theme = useTheme()

    const [mobileNavOpen, setMobileNavOpen] = React.useState(false)
    const isDesktop = useMediaQuery(theme.breakpoints.up(theme.breakpoints.values.desktop))

    React.useEffect(() => {
        if (isDesktop && mobileNavOpen) {
            setMobileNavOpen(false)
        }
    }, [isDesktop])

    const mainNavRef = React.createRef<HTMLElement>()
    const isStuck = useIsStuck(mainNavRef)

    const opacityTransition = theme.transitions.create(
        'opacity',
        {
            duration: theme.transitions.duration.ease,
            easing: theme.transitions.easing.ease,
        }
    )

    const handleClickNavButton = () => {
        setMobileNavOpen(val => !val)
    }

    const [scrollDistanceToMainNav, setScrollDistanceToMainNav] = React.useState(0)

    React.useEffect(() => {
        const mainNavEl = mainNavRef.current

        if (mainNavEl && !isStuck && mobileNavOpen) {
            const setDistanceToMainNav = _.throttle(() => {
                setScrollDistanceToMainNav(
                    Math.max(
                        mainNavEl.getBoundingClientRect().top,
                        0,
                    )
                )
            }, 10)

            document.addEventListener('scroll', setDistanceToMainNav)

            return () => {
                document.removeEventListener('scroll', setDistanceToMainNav)
            }
        }
    }, [isStuck, mobileNavOpen])

    return (
        <React.Fragment>
            {/* UCSF Nav */}
            <Box
                component='nav'
                aria-label='UCSF Main'
                py='10px'
                bgcolor='utilityUCSFNavy.main'
            >
                <Container
                    sx={{
                        display: 'flex',
                        justifyContent: 'space-between',
                    }}
                >
                    <UCSFHomeLink />
                    <Box
                        component='span'
                        sx={(theme) => ({
                            display: 'none',
                            [theme.breakpoints.up('desktop')]: {
                                display: 'inline-block',
                                ml: '29px',
                            }
                        })}
                    >
                        <UCSFNav />
                    </Box>
                </Container>
            </Box>

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
                    height: `${NavBarHeights.condensed}px`,
                    overflow: 'visible',
                    [`& .${NavBarClasses.root}`]: {
                        transition: opacityTransition,
                        position: 'absolute',
                        top: 0,
                        width: '100%',
                    },
                }}
            >
                {!isDesktop && (
                    <MobileNav
                        open={mobileNavOpen}
                        sx={{
                            pt: isStuck ? `${NavBarHeights.default}px` : `${36 + NavBarHeights.default}px`,
                            height: `calc(100vh - ${scrollDistanceToMainNav}px)`,
                        }}
                    />
                )}

                <NavBar
                    variant='default'
                    sx={{
                        opacity: isStuck ? 0 : 1,

                    }}
                    onClickNavButton={handleClickNavButton}
                />

                <NavBar
                    variant='condensed'
                    sx={{
                        opacity: isStuck ? 1 : 0,
                    }}
                    onClickNavButton={handleClickNavButton}
                />


            </Box>
        </React.Fragment>
    )
}