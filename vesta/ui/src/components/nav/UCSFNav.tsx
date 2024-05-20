import * as React from 'react'
import MUILink from '@mui/material/Link';
import Link from 'next/link'
import Box from '@mui/material/Box';


function NavLink({ text, href }: { text: string, href: string }) {
    return (
        <Box
            component='li'
            sx={{
                display: 'inline-block',
                mr: '29px',
            }}
        >
            <MUILink
                href={href}
                tabIndex={0}
                component={Link}
                underline='none'
                color='utilityUCSFLightNavy.main'
            >
                {text}
            </MUILink>
        </Box>
    )
}


export default function UCSFNav() {
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
            <NavLink
                text='About UCSF'
                href='https://www.ucsf.edu/about'
            />
            <NavLink
                text='Search UCSF'
                href='https://www.ucsf.edu/search'
            />
            <NavLink
                text='UCSF Health'
                href='https://www.ucsfhealth.org/'
            />
            <Box
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
                    Donate
                </MUILink>
            </Box>
        </Box>
    )
}