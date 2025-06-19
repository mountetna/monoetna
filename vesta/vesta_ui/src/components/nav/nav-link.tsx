import * as React from 'react'
import Box from '@mui/material/Box';
import { Typography, useTheme } from '@mui/material';

import { TypographyVariant } from '@/lib/utils/types'
import Link from '../link/link';
import { Classes as DLNavClasses } from './dl-nav'

export default function NavLink({
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
            className={DLNavClasses.linkContainer}
            component='li'
        >
            <Link
                className={DLNavClasses.link}
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
