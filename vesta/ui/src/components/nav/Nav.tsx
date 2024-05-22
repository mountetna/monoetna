'use client'

import * as React from 'react'
import Container from '@mui/system/Container'
import ButtonBase from '@mui/material/ButtonBase';
import Box from '@mui/system/Box'
import MUILink from '@mui/material/Link';
import Link from 'next/link'

import UCSFNav from './UCSFNav'
import DLNav from './DLNav'
import UCSFHomeLink from './UCSFHomeLink';
import useIsStuck from '@/lib/utils/useIsStuck';


export default function Nav() {
    const mainNavRef = React.createRef<HTMLElement>()
    const isStuck = useIsStuck(mainNavRef)
    
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
                ref={mainNavRef}
                className={`${isStuck && 'isStuck'}`}
                sx={{
                    position: 'sticky',
                    top: 0,
                    backgroundColor: 'utilityLowlight.main',
                    '&.isStuck': {
                        backgroundColor: 'utilityHighlight.main',

                    },
                }}
            >
                <Container
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
                        sx={(theme) => ({
                            mt: '16px',
                            [theme.breakpoints.up('tablet')]: {
                                mt: '24px',
                            }
                        })}
                    >
                        {/* TODO: logo variants */}
                        LOGO
                    </MUILink>

                    <Box
                        sx={{}}
                    >
                        <DLNav isStuck={isStuck} />

                        <ButtonBase
                            tabIndex={0}
                            sx={(theme) => ({
                                color: 'white',
                                p: '8px',
                                [theme.breakpoints.up('desktop')]: {
                                    display: 'none',
                                }
                            })}
                        >
                            {/* TODO: replace with hamburger icon */}
                            HAMB
                        </ButtonBase>
                    </Box>
                </Container>
            </Box>
        </React.Fragment>
    )
}