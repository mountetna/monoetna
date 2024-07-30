'use client'

import * as React from 'react'
import MUILink from '@mui/material/Link';
import Collapse from '@mui/material/Collapse'
import { TransitionProps } from '@mui/material/transitions';
import Fade from '@mui/material/Fade'
import { SxProps, useTheme } from '@mui/material';
import Box from '@mui/material/Box';
import Link from 'next/link'
import { useRouter } from 'next/navigation';

import LibraryCardButton from '../library-card/library-card-button';
import DLNav, { Classes as DLNavClasses } from './dl-nav';
import UCSFNav from './ucsf-nav';


export default function MobileNav({
    open,
    sx = {},
}: {
    open: boolean,
    sx?: SxProps,
}) {
    const theme = useTheme()

    const animationProps: TransitionProps = {
        in: open,
        easing: theme.transitions.easing.quint,
        timeout: theme.transitions.duration.quint,
    }

    return (
        <Collapse
            component='nav'
            aria-label='Main'
            {...animationProps}
        >
            <Fade
                {...animationProps}
            >
                <Box
                    sx={{
                        bgcolor: 'utilityLowlight.main',
                        px: '8px',
                        py: '8px',
                        ...sx,
                    }}
                >
                    <Box>
                        <DLNav
                            sx={{
                                [`& .${DLNavClasses.link}`]: {
                                    position: 'relative',
                                    color: 'utilityHighlight.main',
                                    left: 0,
                                    '&:hover, &:focus': {
                                        color: 'blue.grade50',
                                        left: '16px',
                                    },
                                    transition: theme.transitions.create(
                                        ['color', 'left'],
                                        {
                                            easing: theme.transitions.easing.quint,
                                            duration: theme.transitions.duration.quint,
                                        },
                                    ),
                                },
                                [`& .${DLNavClasses.libraryCardListItemContainer}`]: {
                                    display: 'none',
                                }
                            }}
                            typography='h2'
                        />
                        <LibraryCardButton />
                    </Box>

                    <UCSFNav />
                </Box>

            </Fade>
        </Collapse>
    )
}