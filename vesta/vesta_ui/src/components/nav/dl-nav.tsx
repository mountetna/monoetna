import * as React from 'react'
import { useRouter } from 'next/navigation';
import Box from '@mui/material/Box';
import { SxProps, Typography, useTheme } from '@mui/material';

import { TypographyVariant } from '@/lib/utils/types'
import { Heights as NavBarHeights } from './nav-bar';
import { useBreakpoint } from '@/lib/utils/responsive';
import Link from '../link/link';
import { toSearchParamsString } from '@/lib/utils/uri';
import { ABOUT_SERACH_PARAMS_KEY, AboutSearchParamsState } from '../about/models';


export enum Classes {
    root = 'data-library-nav',
    linkContainer = 'link-container',
    link = 'link',
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
    const theme = useTheme()

    return (
        <Box
            className={Classes.linkContainer}
            component='li'
        >
            <Link
                className={Classes.link}
                href={href}
                tabIndex={0}
                onClick={onClick}
                sx={{
                    '&:hover, &:focus': {
                        color: 'blue.grade50',
                    },
                    transition: theme.transitions.create(
                        ['color'],
                        {
                            easing: theme.transitions.easing.ease,
                            duration: theme.transitions.duration.ease,
                        },
                    ),
                }}
            >
                <Typography variant={typography}>
                    {text}
                </Typography>
            </Link>
        </Box>
    )
}


export default function DLNav({
    sx = {},
    linkTypography = 'pBody',
    onClickNavLink,
}: {
    sx?: SxProps,
    linkTypography?: TypographyVariant
    onClickNavLink?: () => void,
}) {
    const router = useRouter()

    const breakpoint = useBreakpoint()

    const handleClickNavLink = (event: React.MouseEvent<HTMLAnchorElement>) => {
        event.preventDefault()

        const href = event.currentTarget.href
        const elId = href.split('#')[1]
        const el = document.getElementById(elId)

        router.push(href, { scroll: false })

        if (el) {
            window.scrollTo({
                top: el.offsetTop - NavBarHeights[breakpoint].condensed,
                behavior: 'smooth',
            })
        }

        onClickNavLink && onClickNavLink()
    }

    const aboutSearchParams: AboutSearchParamsState = { index: 0 }
    const aboutHref = '/?' + toSearchParamsString({ [ABOUT_SERACH_PARAMS_KEY]: aboutSearchParams }) + '#about'

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
            <NavLink
                text='About the Library'
                href={aboutHref}
                onClick={handleClickNavLink}
                typography={linkTypography}
            />
            <NavLink
                text='Themes'
                href='/#themes'
                onClick={handleClickNavLink}
                typography={linkTypography}
            />
            <NavLink
                text='Projects'
                href='/#projects'
                onClick={handleClickNavLink}
                typography={linkTypography}
            />
            <NavLink
                text='Contribute Data'
                href='#contribute'
                onClick={handleClickNavLink}
                typography={linkTypography}
            />
            <NavLink
                text='Legacy Version'
                href={'https://datalibraryarchive.ucsf.edu/'}
                onClick={handleClickNavLink}
                typography={linkTypography}
            />
        </Box>
    )
}