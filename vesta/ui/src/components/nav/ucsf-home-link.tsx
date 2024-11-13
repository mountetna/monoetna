import * as React from 'react'
import Box from '@mui/system/Box';

import Link from '../link/link';


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
                <Link
                    href="https://www.ucsf.edu/"
                    tabIndex={0}
                    color='utilityUCSFLightNavy.main'
                >
                    University of California San Francisco
                </Link>
            </Box>
        </Box>
    )
}