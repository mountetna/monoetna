import * as React from 'react'
import MUILink from '@mui/material/Link';
import Link from 'next/link'
import { useRouter } from 'next/navigation';
import Box from '@mui/material/Box';
import { SxProps, Typography } from '@mui/material';

import LibraryCardButton, { Classes as LibraryCardButtonClasses } from '../library-card/library-card-button';
import LibraryCardTray, { Classes as LibraryCardTrayClasses } from '../library-card/library-card-tray';
import { TypographyVariant } from '@/lib/utils/types'
import { Heights as NavBarHeights } from './nav-bar';


export enum Classes {
    root = 'data-library-nav',
    linksContainer = 'links-container',
    linkContainer = 'link-container',
    link = 'link',
    libraryCardListItemContainer = 'library-card-list-item-container',
    libraryCardButton = LibraryCardButtonClasses.root,
    LibraryCardTray = LibraryCardTrayClasses.root,
}


function NavLink({
    text,
    href,
    onClick,
    typography,
}: {
    text: string,
    href: string,
    onClick: (event: React.MouseEvent<HTMLAnchorElement>) => void,
    typography: TypographyVariant,
}) {
    return (
        <Box
            className={Classes.linkContainer}
            component='li'
            sx={{

            }}
        >
            <MUILink
                className={Classes.link}
                href={href}
                tabIndex={0}
                component={Link}
                underline='none'
                onClick={onClick}
            >
                <Typography variant={typography}>
                    {text}
                </Typography>
            </MUILink>
        </Box>
    )
}


export default function DLNav({
    sx = {},
    typography = 'pBody',
}: {
    sx?: SxProps,
    typography?: TypographyVariant
}) {
    const router = useRouter()

    const handleClickNavLink = (event: React.MouseEvent<HTMLAnchorElement>) => {
        event.preventDefault()

        const href = event.currentTarget.href
        const elId = href.split('#')[1]
        const el = document.getElementById(elId)

        router.push(href, { scroll: false })

        if (el) {
            window.scrollTo({
                top: el.offsetTop - NavBarHeights.condensed,
                behavior: 'smooth',
            })
        }
    }

    return (
        <Box
            className={Classes.root}
            component='ol'
            sx={{
                display: 'flex',
                alignItems: 'center',
                p: 0,
                m: 0,
                listStyle: 'none',
                ...sx,
            }}
        >
            <Box
                className={Classes.linksContainer}
            >
                <NavLink
                    text='About the Library'
                    href='#about'
                    onClick={handleClickNavLink}
                    typography={typography}
                />
                <NavLink
                    text='Themes'
                    href='#themes'
                    onClick={handleClickNavLink}
                    typography={typography}
                />
                <NavLink
                    text='Projects'
                    href='#projects'
                    onClick={handleClickNavLink}
                    typography={typography}
                />
                {/* <NavLink
                    text='Contibute'
                    href='#'
                    onClick={handleClickNavLink}
                    typography={typography}
                /> */}
            </Box>
            <Box
                className={Classes.libraryCardListItemContainer}
                component='li'
            >
                <LibraryCardButton />
                <LibraryCardTray />
            </Box>
        </Box>
    )
}