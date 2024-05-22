import * as React from 'react'
import MUILink from '@mui/material/Link';
import Link from 'next/link'
import Box from '@mui/material/Box';
import LibraryCardButton from '../library-card/LibraryCardButton';
import LibraryCardTray from '../library-card/LibraryCardTray';


function NavLink({ text, href, isStuck }: { text: string, href: string, isStuck: boolean }) {
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
                color={isStuck ? 'utilityLowlight.main' : 'utilityHighlight.main'}
                typography='pBodyMediumWt'
            >
                {text}
            </MUILink>
        </Box>
    )
}


export default function DLNav({ isStuck }: { isStuck: boolean }) {
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
                    isStuck={isStuck}
                />
                <NavLink
                    text='Themes'
                    href='/themes'
                    isStuck={isStuck}
                />
                <NavLink
                    text='Projects'
                    href='/projects'
                    isStuck={isStuck}
                />
                {/* <NavLink
                    text='Contibute'
                    href='#'
                    isStuck={isStuck}
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