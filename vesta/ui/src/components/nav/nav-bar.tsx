'use client'

import * as React from 'react'
import Container from '@mui/system/Container'
import ButtonBase from '@mui/material/ButtonBase';
import Box from '@mui/system/Box'
import MUILink from '@mui/material/Link';
import { SxProps, useTheme } from '@mui/material';
import Link from 'next/link'
import Image from 'next/image';

import DLNav from './dl-nav'

import hamburgerIconLightSrc from '/public/images/icons/hamburger-menu-icon-light.svg'
import hamburgerIconDarkSrc from '/public/images/icons/hamburger-menu-icon-dark.svg'
import logoLightSrc from '/public/images/logo/logo-light.svg'
import logoDarkSrc from '/public/images/logo/logo-dark.svg'
import logoWordmarkBottomLightSrc from '/public/images/logo/logo-wordmark-bottom-light.svg'


export enum Classes {
    root = 'dl-navbar',
    mobileLogoContainer = 'logo-container-mobile',
    logoBase = 'logo',
    logoNoWordmarkDark = 'logo-no-wordmark-dark',
    logoNoWordmarkLight = 'logo-no-wordmark-light',
    logoWordmarkBottomLight = 'logo-wordmark-bottom-light',
}

export enum Heights {
    default = 84.9,
    condensed = 64,
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

    return (
        <Box
            className={Classes.root}
            sx={{
                px: '8px',
                bgcolor: variant === 'condensed' ? 'utilityHighlight.main' : 'transparent',
                height: variant === 'condensed' ? Heights.condensed : Heights.default,
                ...sx,
            }}
        >
            <Container
                sx={{
                    display: 'flex',
                    justifyContent: 'space-between',
                    alignItems: 'center',
                    height: '100%',
                    [theme.breakpoints.up('tablet')]: {
                        // alignItems: 'flex-start',
                        // pl: '32px',
                        // pr: '16px',
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
                        mt: variant === 'condensed' ? 'unset' : '16px',
                        [`& .${Classes.logoBase}`]: {
                            width: 'auto',
                            height: '100%',
                            '& img': {
                                width: 'auto',
                                height: '100%',
                            },
                        },
                        [theme.breakpoints.up('tablet')]: {
                            mt: '24px',
                            width: '101px',
                            height: '154px',
                        }
                    })}
                >
                    <Box
                        className={Classes.mobileLogoContainer}
                        sx={{
                            display: 'flex',
                            position: 'relative',
                            width: variant === 'condensed' ? '53.1px' : '74.09px',
                            height: variant === 'condensed' ? '53.1px' : '74.09px',
                            '& > *': {
                                display: 'flex',
                            },
                            [theme.breakpoints.up('tablet')]: {
                                display: 'none',
                            }
                        }}
                    >
                        <Box
                            className={[Classes.logoBase, Classes.logoNoWordmarkLight].join(' ')}
                            sx={{
                                opacity: variant === 'condensed' ? 0 : 1,
                            }}
                        >
                            <Image
                                src={logoLightSrc}
                                alt='Data Library Logo'
                            />
                        </Box>

                        <Box
                            className={[Classes.logoBase, Classes.logoNoWordmarkDark].join(' ')}
                            sx={{
                                position: 'absolute',
                                top: 0,
                                left: 0,
                                opacity: variant === 'condensed' ? 1 : 0,
                            }}
                        >
                            <Image
                                src={logoDarkSrc}
                                alt='Data Library Logo'
                            />
                        </Box>
                    </Box>

                    {/* <Box
                            className='logo logo-wordmark'
                            sx={{
                                opacity: '0',
                                [theme.breakpoints.up('tablet')]: {
                                    opacity: '1',
                                }
                            }}
                        >
                            <Image
                                src={logoWordmarkBottomLightSrc}
                                alt='Data Library Logo with wordmark'
                            />
                        </Box> */}
                </MUILink>

                <Box>
                    <DLNav
                        sx={{
                            display: 'none',
                        }}
                    />

                    <ButtonBase
                        tabIndex={0}
                        onClick={() => onClickNavButton()}
                        sx={(theme) => ({
                            color: 'utilityHighlight.main',
                            width: '48px',
                            height: '48px',
                            mx: '9px',
                            my: '8px',
                            position: 'relative',
                            '& img': {
                                position: 'absolute',
                            },
                            '& .icon-dark': {
                                opacity: variant === 'condensed' ? 1 : 0,
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
        </Box>
    )
}