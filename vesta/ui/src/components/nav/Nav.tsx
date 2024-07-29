'use client'

import * as React from 'react'
import Container from '@mui/system/Container'
import ButtonBase from '@mui/material/ButtonBase';
import Box from '@mui/system/Box'
import MUILink from '@mui/material/Link';
import Link from 'next/link'
import Image from 'next/image';

import UCSFNav from './ucsf-nav'
import DLNav from './dl-nav'
import UCSFHomeLink from './ucsf-home-link';
import useIsStuck from '@/lib/utils/css';
import { useTheme } from '@mui/material/styles';

import hamburgerIconLightSrc from '/public/images/icons/hamburger-menu-icon-light.svg'
import hamburgerIconDarkSrc from '/public/images/icons/hamburger-menu-icon-dark.svg'
import logoLightSrc from '/public/images/logo/logo-light.svg'
import logoDarkSrc from '/public/images/logo/logo-dark.svg'
import logoWordmarkBottomLightSrc from '/public/images/logo/logo-wordmark-bottom-light.svg'
import MobileNav from './mobile-nav';


export default function Nav() {
    const [mobileNavOpen, setMobileNavOpen] = React.useState(false)

    const mainNavRef = React.createRef<HTMLElement>()
    // TODO: fix glitchiness
    // const isStuck = useIsStuck(mainNavRef)
    const isStuck = true
    const theme = useTheme()
    const transition = theme.transitions.create(
        // ['opacity', 'color', 'background-color'],
        ['all'],
        {
            duration: theme.transitions.duration.quint,
            easing: theme.transitions.easing.quint,
        }
    )

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
                component='nav'
                aria-label='Main'
                className={`${isStuck ? 'stuck' : ''}`}
                ref={mainNavRef}
                zIndex='appBar'
                sx={(theme) => ({
                    position: 'sticky',
                    top: 0,
                    width: '100%',
                    backgroundColor: 'utilityLowlight.main',
                    transition: transition,
                    px: '8px',
                    '&.stuck': {
                        // position: 'fixed',
                        top: 0,
                        backgroundColor: 'utilityHighlight.main',
                        '.logo-no-wordmark-dark, .icon-dark': {
                            opacity: '1',
                        },
                        '.logo-wordmark, .logo-no-wordmark-light, .icon-light': {
                            opacity: '0',
                        },
                        '.home-link': {
                            mt: 0,
                            width: '49px',
                            height: '49px',
                            [theme.breakpoints.up('tablet')]: {
                                width: '77px',
                                height: '77px',
                            }
                        },
                        '.dl-nav-container': {
                            [theme.breakpoints.up('tablet')]: {
                                alignItems: 'center',
                                py: '16px',
                            },
                        },
                    },

                })}
            >
                <Container
                    className='dl-nav-container'
                    sx={{
                        display: 'flex',
                        justifyContent: 'space-between',
                        alignItems: 'center',
                        // transition: transition,
                        [theme.breakpoints.up('tablet')]: {
                            alignItems: 'flex-start',
                            pl: '32px',
                            pr: '16px',
                        },
                        [theme.breakpoints.up('desktop')]: {
                        },
                    }}
                >
                    <MUILink
                        href="/"
                        tabIndex={0}
                        component={Link}
                        underline='none'
                        className='home-link'
                        sx={(theme) => ({
                            position: 'relative',
                            width: '77px',
                            height: '77px',
                            '& .logo': {
                                transition: transition,
                                position: 'absolute',
                                width: 'auto',
                                height: '100%',
                                'img': {
                                    width: 'auto',
                                    height: '100%',
                                },
                            },
                            color: 'utilityHighlight.main',
                            mt: '16px',
                            [theme.breakpoints.up('tablet')]: {
                                mt: '24px',
                                width: '101px',
                                height: '154px',
                            }
                        })}
                    >
                        <Box
                            className='logo logo-no-wordmark-light'
                            sx={(theme) => ({
                                [theme.breakpoints.up('tablet')]: {
                                    opacity: '0',
                                }
                            })}
                        >
                            <Image
                                src={logoLightSrc}
                                alt='Data Library Logo'
                            />
                        </Box>
                        <Box
                            className='logo logo-no-wordmark-dark'
                            sx={(theme) => ({
                                opacity: '0',
                            })}
                        >
                            <Image
                                src={logoDarkSrc}
                                alt='Data Library Logo'
                            />
                        </Box>
                        <Box
                            className='logo logo-wordmark'
                            sx={(theme) => ({
                                opacity: '0',
                                [theme.breakpoints.up('tablet')]: {
                                    opacity: '1',
                                }
                            })}
                        >
                            <Image
                                src={logoWordmarkBottomLightSrc}
                                alt='Data Library Logo with wordmark'
                            />
                        </Box>
                    </MUILink>

                    <Box>
                        <DLNav isStuck={isStuck} />

                        <ButtonBase
                            tabIndex={0}
                            onClick={() => setMobileNavOpen(!mobileNavOpen)}
                            sx={(theme) => ({
                                color: 'utilityHighlight.main',
                                width: '48px',
                                height: '48px',
                                mx: '9px',
                                my: '8px',
                                position: 'relative',
                                '& img': {
                                    position: 'absolute',
                                    transition: transition,
                                },
                                '& .icon-dark': {
                                    opacity: '0',
                                },
                                [theme.breakpoints.up('tablet')]: {
                                    m: '16px',
                                },
                                [theme.breakpoints.up('desktop')]: {
                                    display: 'none',
                                },
                            })}
                        >
                            <Image
                                src={hamburgerIconLightSrc}
                                alt='Hamburger Nav Icon'
                                className='hamburger-nav-icon icon-light'
                            />
                            <Image
                                src={hamburgerIconDarkSrc}
                                alt='Hamburger Nav Icon'
                                className='hamburger-nav-icon icon-dark'
                            />
                        </ButtonBase>
                    </Box>
                </Container>

                <MobileNav
                    open={mobileNavOpen}
                />
            </Box>
        </React.Fragment>
    )
}