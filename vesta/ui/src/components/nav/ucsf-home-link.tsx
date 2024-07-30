import * as React from 'react'
import MUILink from '@mui/material/Link';
import Link from 'next/link'
import Box from '@mui/system/Box';


export default function UCSFHomeLink() {
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
                component='li'
                sx={{

                }}
            >
                <MUILink
                    href="https://www.ucsf.edu/"
                    tabIndex={0}
                    component={Link}
                    underline='none'
                    color='utilityUCSFLightNavy.main'
                >
                    University of California San Francisco
                </MUILink>
            </Box>
        </Box>
    )
}