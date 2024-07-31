'use client'

import * as React from 'react';
import Container from '@mui/system/Container'
import Box from '@mui/system/Box'
import MUILink from '@mui/material/Link';
import Link from 'next/link'
import Image from 'next/image';

import logoWordmarkRightLightSrc from '/public/images/logo/logo-wordmark-right-light.svg'
import Copyright from './legal/copyright';


function NavLink({ text, href }: { text: string, href: string }) {
    return (
        <Box
            component='li'
            sx={{ display: 'block' }}
        >
            <MUILink
                href={href}
                tabIndex={0}
                component={Link}
                underline='none'
                color='ground.grade100'
                sx={(theme) => ({
                    typography: 'pBody',
                    [theme.breakpoints.up('tablet')]: {
                        typography: 'pMedium',
                    }
                })}
            >
                {text}
            </MUILink>
        </Box>
    )
}

export default function Footer() {
    return (
        <Container
            component='footer'
        >
            <Box
                sx={(theme) => ({
                    p: '16px',
                    borderRadius: '30px',
                    color: 'ground.grade100',
                    bgcolor: 'ground.grade10',
                    [theme.breakpoints.up('tablet')]: {
                        px: '44px',
                        py: '16px',
                    },
                    [theme.breakpoints.up('desktop')]: {
                        px: '24px',
                        py: '16px',
                    },
                })}
            >
                <Box
                    sx={(theme) => ({
                        display: 'grid',
                        gridTemplateColumns: 'repeat(6, 1fr)',
                        gridTemplateRows: 'auto',
                        mb: '8rem', // TODO: is this right? should the whole footer be deterministically-sized?
                        [theme.breakpoints.up('tablet')]: {
                            gridTemplateColumns: 'repeat(12, 1fr)',
                            mb: '12rem',
                        },
                    })}
                >
                    {/* TODO: replace */}
                    <Box
                        sx={(theme) => ({
                            color: 'utilityHighlight.main',
                            '& svg': {
                                width: 'auto',
                                height: '100%',
                            },
                            gridColumn: 'span 6',
                            mb: '16px',
                            [theme.breakpoints.up('tablet')]: {
                                mb: '0',
                            },
                            [theme.breakpoints.up('desktop')]: {
                                gridColumn: 'span 4',
                            },
                        })}
                    >
                        {/* TODO: switch to right aligned wordmark */}
                        <Image
                            src={logoWordmarkRightLightSrc}
                            alt='adsa'
                        />
                    </Box>
                    <Box
                        component='nav'
                        aria-label='UCSF Legal'
                        sx={(theme) => ({
                            gridColumn: 'span 3',
                            [theme.breakpoints.up('desktop')]: {
                                gridColumn: 'span 2',
                            },
                        })}
                    >
                        <Box
                            component='ol'
                            sx={(theme) => ({
                                display: 'block',
                                p: 0,
                                m: 0,
                                listStyle: 'none',
                                '& > *:not(:last-child)': {
                                    mb: '10px',
                                },
                                [theme.breakpoints.up('tablet')]: {
                                    py: '16px',
                                },
                            })}
                        >
                            <Box
                                component='li'
                                sx={(theme) => ({
                                    display: 'block',
                                    typography: 'pBodyBoldWt',
                                    [theme.breakpoints.up('tablet')]: {
                                        typography: 'pMediumBoldWt',
                                    }
                                })}
                            >
                                General
                            </Box>
                            <NavLink
                                text='Accessibility'
                                href='https://www.ucsf.edu/accessibility-resources'
                            />
                            <NavLink
                                text='Privacy Policy'
                                href='https://www.ucsf.edu/website-privacy-policy'
                            />
                            <NavLink
                                text='Terms of Use'
                                href='https://websites.ucsf.edu/website-terms-use'
                            />
                            <NavLink
                                text='A-Z Website List'
                                href='https://websites.ucsf.edu/azlist'
                            />
                        </Box>
                    </Box>
                </Box>
                
                <Copyright />
            </Box>
        </Container>
    )
}