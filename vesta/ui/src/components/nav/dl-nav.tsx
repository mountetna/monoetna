import * as React from 'react'
import MUILink from '@mui/material/Link';
import Link from 'next/link'
import Box from '@mui/material/Box';

import LibraryCardButton from '../library-card/library-card-button';
import LibraryCardTray from '../library-card/library-card-tray';
import { useRouter } from 'next/navigation';


function NavLink({ text, href, isStuck, onClick }: {
    text: string,
    href: string,
    isStuck: boolean,
    onClick: (event: React.MouseEvent<HTMLAnchorElement>) => void,
}) {
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
                onClick={onClick}
            >
                {text}
            </MUILink>
        </Box>
    )
}


export default function DLNav({ isStuck }: { isStuck: boolean }) {
    const router = useRouter()

    const handleClickNavLink = (event: React.MouseEvent<HTMLAnchorElement>) => {
        event.preventDefault()

        const href = event.currentTarget.href
        const elId = href.split('#')[1]
        const el = document.getElementById(elId)

        router.push(href, { scroll: false })
        window.scrollTo({
            top: el?.offsetTop,
            behavior: 'smooth',
        })
    }

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
                        [theme.breakpoints.up('desktop')]: {
                            mr: '54px',
                        },
                    },
                    [theme.breakpoints.up('desktop')]: {
                        display: 'inline-block',
                        pl: '30px',
                        py: '18px'
                    },
                })}
            >
                <NavLink
                    text='About the Library'
                    href='#about'
                    isStuck={isStuck}
                    onClick={handleClickNavLink}
                />
                <NavLink
                    text='Themes'
                    href='#'
                    isStuck={isStuck}
                    onClick={handleClickNavLink}
                />
                <NavLink
                    text='Projects'
                    href='#'
                    isStuck={isStuck}
                    onClick={handleClickNavLink}
                />
                <NavLink
                    text='Contibute'
                    href='#'
                    isStuck={isStuck}
                    onClick={handleClickNavLink}
                />
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