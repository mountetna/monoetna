import * as React from 'react'
import { SxProps, Typography } from '@mui/material';

import { TypographyVariant } from '@/lib/utils/types';


export default function Copyright({
    typography = 'pXS',
    sx = {},
}: {
    typography?: TypographyVariant,
    sx?: SxProps,
}) {
    return (
        <Typography
            variant={typography}
            sx={{
                ...sx
            }}
        >
            © {(new Date()).getFullYear()} The Regents of the University of California
            <br />
            All rights reserved.
        </Typography>
    )
}