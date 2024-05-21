import * as React from 'react'
import MUILink from '@mui/material/Link';
import Link from 'next/link'
import Box from '@mui/material/Box';
import LibraryCardButton from '../library-card/LibraryCardButton';
import LibraryCardTray from '../library-card/LibraryCardTray';


function NavLink({ text, href }: { text: string, href: string }) {
    return (
        <Box
            component='li'
            sx={{
                display: 'inline-block',
                py: '8px',
            }}
        >
            <MUILink
                href={href}
                tabIndex={0}
                component={Link}
                underline='none'
                color='utilityHighlight.main'
                typography='pBodyMediumWt'
            >
                {text}
            </MUILink>
        </Box>
    )
}


export default function DLNav() {
    return (
        <Box
            component='ol'
            sx={{
                display: 'inline-flex',
                alignItems: 'center',
                p: 0,
                m: 0,
                listStyle: 'none',
                typography: 'pBody',
            }}
        >
            <Box
                sx={(theme) => ({
                    display: 'none',
                    '& > *': {
                        mr: '29px',
                    },
                    [theme.breakpoints.up('desktop')]: {
                        display: 'inline-block',
                    },
                })}
            >
                <NavLink
                    text='About the Library'
                    href='/about'
                />
                <NavLink
                    text='Themes'
                    href='/themes'
                />
                <NavLink
                    text='Projects'
                    href='/projects'
                />
                {/* <NavLink
                    text='Contibute'
                    href='#'
                /> */}
            </Box>
            <Box
                component='li'
                sx={(theme) => ({
                    display: 'none',
                    [theme.breakpoints.up('tablet')]: {
                        display: 'inline-block',
                    },
                })}
            >
                <Box
                    sx={(theme) => ({
                        display: 'none',
                        [theme.breakpoints.up('tablet')]: {
                            display: 'inline-block',
                        },
                        [theme.breakpoints.up('desktop')]: {
                            display: 'none',
                        },
                    })}
                >
                    <LibraryCardButton />
                </Box>
                <Box
                    sx={(theme) => ({
                        display: 'none',
                        [theme.breakpoints.up('desktop')]: {
                            display: 'inline-block',
                            mr: '16px',
                        },
                    })}
                >
                    <LibraryCardTray />
                </Box>
            </Box>
        </Box>
    )
}