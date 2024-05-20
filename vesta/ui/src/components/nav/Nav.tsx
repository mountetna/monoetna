'use client'
// needed to use sx with function param
// TODO: figure out why necessary

import * as React from 'react'
import Container from '@mui/system/Container'
import ButtonBase from '@mui/material/ButtonBase';
import Box from '@mui/system/Box'
import MUILink from '@mui/material/Link';
import Link from 'next/link'

import UCSFNav from './UCSFNav'
import DLNav from './DLNav'
import UCSFHomeLink from './UCSFHomeLink';


export default function Nav() {
    return (
        <nav>
            <Box
                sx={{
                    py: '10px',
                    backgroundColor: 'utilityUCSFNavy.main',
                }}
            >
                <Container
                    maxWidth="desktopLg"
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
                            [theme.breakpoints.up('tablet')]: {
                                display: 'inline-block',
                                ml: '29px',
                            }
                        })}
                    >
                        <UCSFNav />
                    </Box>
                </Container>
            </Box>
            <Box sx={{ background: 'black' }}>
                <Container
                    maxWidth="desktopLg"
                    sx={{
                        display: 'flex',
                        justifyContent: 'space-between',
                        alignItems: 'flex-start',
                    }}
                >
                    <MUILink
                        href="/"
                        tabIndex={0}
                        component={Link}
                        underline='none'
                    >
                        {/* TODO: logo variants */}
                        LOGO
                    </MUILink>

                    <ButtonBase
                        // variant='text'
                        tabIndex={0}
                        sx={(theme) => ({
                            [theme.breakpoints.up('desktop')]: {
                                display: 'none',
                            }
                        })}
                    >
                        {/* TODO: replace with hamburger icon */}
                        HAMB
                    </ButtonBase>
                </Container>
                <Box>

                    <DLNav />
                </Box>
            </Box>
        </nav>
    )
}