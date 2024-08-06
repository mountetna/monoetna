'use client'

import * as React from 'react'
import Container from '@mui/system/Container'
import ButtonBase from '@mui/material/ButtonBase';
import Box from '@mui/system/Box'
import MUILink from '@mui/material/Link';
import { Breakpoint, Collapse, SxProps, useTheme } from '@mui/material';
import Link from 'next/link'
import Image from 'next/image';

import DLNav, { Classes as DLNavClasses } from './dl-nav'
import { useBreakpoint } from '@/lib/utils/responsive';
import LibraryCardButton, { Classes as LibraryCardButtonClasses } from '../library-card/library-card-button';
import { LibraryCardModal } from '../library-card/library-card-modal';
import { useUser } from '../user/context';
import LibraryCardTray, { Classes as LibraryCardTrayClasses } from '../library-card/library-card-tray';

import hamburgerIconLightSrc from '/public/images/icons/hamburger-menu-icon-light.svg'
import hamburgerIconDarkSrc from '/public/images/icons/hamburger-menu-icon-dark.svg'
import logoLightSrc from '/public/images/logo/logo-light.svg'
import logoDarkSrc from '/public/images/logo/logo-dark.svg'
import wordmarkBottomLightSrc from '/public/images/logo/wordmark-bottom-light.svg'


export enum Classes {
    root = 'dl-navbar',
    logoContainer = 'logo-container',
    logo = 'logo',
    logoWordmark = 'logo-wordmark',
}

interface BreakpointHeights {
    default: number
    condensed: number
}

export const Heights: Record<Breakpoint, BreakpointHeights> = {
    mobile: { default: 84.9, condensed: 64 },
    tablet: { default: 174, condensed: 95 },
    desktop: { default: 173.77, condensed: 92.95 },
    desktopLg: { default: 173.77, condensed: 92.95 },
}


export default function NavBar({
    variant,
    onClickNavButton,
    sx = {},
}: {
    variant: 'default' | 'condensed',
    onClickNavButton: () => void,
    sx?: SxProps
}) {
    const theme = useTheme()
    const transition = theme.transitions.create(
        'all',
        {
            duration: theme.transitions.duration.ease,
            easing: theme.transitions.easing.ease,
        },
    )

    const breakpoint = useBreakpoint()
    const isTablet = breakpoint === 'tablet'
    const isDesktop = ['desktop', 'desktopLg'].includes(breakpoint)

    const user = useUser()

    const [libraryCardModalOpen, setLibraryCardModalOpen] = React.useState(false)

    const [libraryCardTrayOpen, setLibraryCardTrayOpen] = React.useState(false)
    const libraryCardTrayRef = React.useRef<HTMLElement>()
    const libraryCardTrayWidth = libraryCardTrayRef.current?.offsetWidth

    // Close Library Card Tray when scroll all the way up
    React.useEffect(() => {
        if (variant === 'default') {
            setLibraryCardTrayOpen(false)
        }
    }, [variant])

    return (
        <Collapse
            in={variant === 'default'}
            collapsedSize={Heights[breakpoint].condensed}
            easing={theme.transitions.easing.ease}
            timeout={theme.transitions.duration.ease}
            className={Classes.root}
            sx={{
                bgcolor: variant === 'condensed' ? 'utilityHighlight.main' : 'transparent',
                overflow: 'visible',
                display: 'flex',
                transition,
                ...sx,
            }}
        >
            <Container
                sx={{
                    display: 'flex',
                    justifyContent: 'space-between',
                    alignItems: 'center',
                    position: 'relative',
                    transition,
                    [theme.breakpoints.up('tablet')]: {
                        pl: '32px',
                        pr: '16px',
                        py: variant === 'condensed' ? '10px' : '0px',
                    },
                    [theme.breakpoints.up('desktop')]: {
                        py: variant === 'condensed' ? '16px' : '0px',
                    },
                }}
            >
                <MUILink
                    href="/"
                    tabIndex={0}
                    component={Link}
                    underline='none'
                    className='home-link'
                    sx={{
                        mt: variant === 'condensed' ? '0px' : '16px',
                        [`& .${Classes.logo}`]: {
                            width: 'auto',
                            height: '100%',
                            '& > img': {
                                width: 'auto',
                                height: '100%',
                            },
                        },
                        transition,
                    }}
                >
                    <Box
                        className={Classes.logoContainer}
                        sx={{
                            display: 'flex',
                            position: 'relative',
                            width: variant === 'condensed' ? '53.1px' : '74.09px',
                            height: variant === 'condensed' ? '53.1px' : '74.09px',
                            transition,
                            '& > *': {
                                display: 'flex',
                            },
                            [theme.breakpoints.up('tablet')]: {
                                width: variant === 'condensed' ? '75px' : '99px',
                                height: variant === 'condensed' ? '75px' : '99px',
                            },
                            [theme.breakpoints.up('desktop')]: {
                                width: variant === 'condensed' ? '60.95px' : '99px',
                                height: variant === 'condensed' ? '60.95px' : '99px',
                            },
                        }}
                    >
                        <Box
                            className={Classes.logo}
                            sx={{
                                opacity: variant === 'condensed' ? 0 : 1,
                                transition,
                                position: 'relative',
                            }}
                        >
                            <Image
                                src={logoLightSrc}
                                alt='Data Library Logo'
                            />

                            <Box
                                sx={{
                                    display: 'none',
                                    [theme.breakpoints.up('tablet')]: {
                                        position: 'absolute',
                                        top: '100%',
                                        display: 'flex',
                                        width: '100%',
                                        pt: '13px',
                                        '& > img': {
                                            width: '100%',
                                            height: 'auto',
                                        }
                                    },
                                }}
                            >
                                <Image
                                    className={Classes.logoWordmark}
                                    src={wordmarkBottomLightSrc}
                                    alt='Data Library Logo Wordmark'
                                />
                            </Box>
                        </Box>

                        <Box
                            className={Classes.logo}
                            sx={{
                                position: 'absolute',
                                top: 0,
                                left: 0,
                                opacity: variant === 'condensed' ? 1 : 0,
                                transition,
                            }}
                        >
                            <Image
                                src={logoDarkSrc}
                                alt='Data Library Logo'
                            />
                        </Box>
                    </Box>
                </MUILink>

                <Box
                    sx={{
                        display: 'flex',
                        px: '9px',
                        py: '8px',
                        gap: variant === 'condensed' ? '16px' : '34px',
                        alignItems: 'center',
                        transition,
                        [theme.breakpoints.up('desktop')]: {
                            pl: '30px',
                            pr: '16px',
                            py: '0px',
                            gap: '54px',
                        },
                    }}
                >
                    <DLNav
                        sx={{
                            display: 'none',
                            [theme.breakpoints.up('desktop')]: {
                                display: 'flex',
                                gap: '54px',
                                [`.${DLNavClasses.link}`]: {
                                    color: variant === 'condensed' ? 'ground.grade10' : 'utilityHighlight.main',
                                    '&:hover, &:focus': {
                                        color: 'blue.grade50',
                                    },
                                },
                            },
                        }}
                    />

                    {isTablet && <LibraryCardButton
                        isLoggedIn={user !== null}
                        onClick={() => setLibraryCardModalOpen(val => !val)}
                    />}

                    {isTablet && user && (
                        <LibraryCardModal
                            open={libraryCardModalOpen}
                            handleSetOpen={setLibraryCardModalOpen}
                            user={user}
                        />
                    )}

                    {isDesktop && (
                        <Box
                            sx={{
                                width: `${libraryCardTrayWidth}px`,
                                height: '1px',
                                '& > *': {
                                    position: 'absolute',
                                    top: '0px',
                                    opacity: libraryCardTrayWidth === undefined ? 0 : 1,
                                    gap: variant === 'condensed' ? '15px' : '26px',
                                },
                            }}
                        >
                            <LibraryCardTray
                                ref={libraryCardTrayRef}
                                open={libraryCardTrayOpen}
                                onSetOpen={setLibraryCardTrayOpen}
                                user={user}
                            />
                        </Box>
                    )}

                    <ButtonBase
                        tabIndex={0}
                        onClick={() => onClickNavButton()}
                        sx={(theme) => ({
                            color: 'utilityHighlight.main',
                            width: '48px',
                            height: '48px',
                            position: 'relative',
                            transition,
                            '& img': {
                                position: 'absolute',
                            },
                            '& .icon-dark': {
                                opacity: variant === 'condensed' ? 1 : 0,
                                transition,
                            },
                            [theme.breakpoints.up('tablet')]: {
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
        </Collapse >
    )
}