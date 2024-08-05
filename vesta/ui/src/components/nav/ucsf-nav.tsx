import * as React from 'react'
import MUILink from '@mui/material/Link';
import Link from 'next/link'
import Box from '@mui/material/Box';
import { SxProps, Typography } from '@mui/material';

import { TypographyVariant } from '@/lib/utils/types';
import Copyright from '../legal/copyright';


export enum Classes {
    root = 'ucsf-nav',
    linkContainer = 'link-container',
    link = 'link',
    donateLinkContainer = 'donate-link-container',
}


function NavLink({
    text,
    href,
    typography,
}: {
    text: string,
    href: string,
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
            >
                <Typography variant={typography}>
                    {text}
                </Typography>
            </MUILink>
        </Box>
    )
}


export default function UCSFNav({
    sx = {},
    linkTypography = 'pBody',
    buttonTypography = 'pBodyMediumWt',
}: {
    sx?: SxProps,
    linkTypography?: TypographyVariant,
    buttonTypography?: TypographyVariant,
}) {
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
                text='About UCSF'
                href='https://www.ucsf.edu/about'
                typography={linkTypography}
            />
            <NavLink
                text='Search UCSF'
                href='https://www.ucsf.edu/search'
                typography={linkTypography}
            />
            <NavLink
                text='UCSF Health'
                href='https://www.ucsfhealth.org/'
                typography={linkTypography}
            />
            <Box
                className={Classes.donateLinkContainer}
                component='li'
                sx={{
                    display: 'inline-block'
                }}
            >
                <MUILink
                    href='https://giving.ucsf.edu/'
                    tabIndex={0}
                    component={Link}
                    underline='none'
                    sx={{
                        display: 'inline-block',
                        px: '16px',
                        py: '8px',
                        borderRadius: '40px',
                        color: 'utilityUCSFNavy.main',
                        backgroundColor: 'utilityUCSFLightNavy.main',
                    }}
                >
                    <Typography variant={buttonTypography}>
                        Donate
                    </Typography>
                </MUILink>
            </Box>
        </Box>
    )
}