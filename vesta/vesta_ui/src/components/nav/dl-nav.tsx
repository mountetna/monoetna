import * as React from 'react'
import { useRouter } from 'next/navigation';
import Box from '@mui/material/Box';
import { SxProps, Typography, useTheme } from '@mui/material';

import { Heights as NavBarHeights } from './nav-bar';
import { useBreakpoint } from '@/lib/utils/responsive';
import RelatedResourcesMenu from './related-resources-menu';
import NavLink from './nav-link';
import { toSearchParamsString } from '@/lib/utils/uri';
import { ABOUT_SERACH_PARAMS_KEY, AboutSearchParamsState } from '../about/models';
import { TypographyVariant } from '@/lib/utils/types'


export enum Classes {
    root = 'data-library-nav',
    linkContainer = 'link-container',
    link = 'link',
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

        if (href.startsWith('/')) router.push(href, { scroll: false })
        else window.open(href, '_blank');

        if (el) {
            window.scrollTo({
                top: el.offsetTop - NavBarHeights[breakpoint].condensed,
                behavior: 'smooth',
            })
        }

        onClickNavLink && onClickNavLink()
    }

    const [anchorEl, setAnchorEl] = React.useState<HTMLAnchorElement|null>(null);
    const handleClick = (event: React.MouseEvent<HTMLAnchorElement>) => setAnchorEl(event.currentTarget);
    const handleClose = () => setAnchorEl(null);

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
                text='Related Resources'
                href={'/'}
                onClick={handleClick}
                typography={linkTypography}
            />
            <RelatedResourcesMenu
              anchorEl={anchorEl}
              typography={linkTypography}
              onClick={handleClickNavLink}
              onClose={handleClose}
            />
        </Box>
    )
}
