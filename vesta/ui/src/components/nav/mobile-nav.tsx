'use client'

import * as React from 'react'
import Collapse from '@mui/material/Collapse'
import { TransitionProps } from '@mui/material/transitions';
import Fade from '@mui/material/Fade'
import { SxProps, useTheme } from '@mui/material';
import Box from '@mui/material/Box';

import LibraryCardButton, { Classes as LibraryCardButtonClasses } from '../library-card/library-card-button';
import DLNav, { Classes as DLNavClasses } from './dl-nav';
import UCSFNav, { Classes as UCSFNavClasses } from './ucsf-nav';
import Copyright from '../legal/copyright';
import { LibraryCardModal } from '../library-card/library-card-modal';
import { useUser } from '../user/context';


export default function MobileNav({
    open,
    sx = {},
}: {
    open: boolean,
    sx?: SxProps,
}) {
    const theme = useTheme()
    const user = useUser()

    const [libraryCardModalOpen, setLibraryCardModalOpen] = React.useState(false)

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
                        px: '16px',
                        py: '16px',
                        display: 'flex',
                        flexDirection: 'column',
                        justifyContent: 'space-between',
                        gap: '32px',
                        pb: '106px',
                        ...sx,
                    }}
                >

                    <Box
                        sx={{
                            display: 'flex',
                            flexDirection: 'column',
                            gap: '24px',
                            [`& .${LibraryCardButtonClasses.root}`]: {
                                width: 'fit-content',
                            },
                        }}
                    >
                        <DLNav
                            sx={{
                                [`.${DLNavClasses.linksContainer}`]: {
                                    display: 'flex',
                                    flexDirection: 'column',
                                    alignItems: 'start',
                                    gap: '16px',
                                },
                                [`.${DLNavClasses.link}`]: {
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
                                [`.${DLNavClasses.libraryCardListItemContainer}`]: {
                                    display: 'none',
                                }
                            }}
                            typography='h2'
                        />

                        <LibraryCardButton
                            onClick={() => setLibraryCardModalOpen(val => !val)}
                        />

                        {user && (
                            <LibraryCardModal
                                open={libraryCardModalOpen}
                                handleSetOpen={setLibraryCardModalOpen}
                                user={user}
                            />
                        )}
                    </Box>

                    <Box
                        sx={{
                            display: 'flex',
                            flexDirection: 'column',
                            gap: '16px',
                            color: 'utilityHighlight.main',
                        }}
                    >
                        <UCSFNav
                            sx={{
                                display: 'flex',
                                flexDirection: 'column',
                                alignItems: 'start',
                                gap: '16px',
                                [`.${UCSFNavClasses.link}`]: {
                                    color: 'utilityHighlight.main',
                                },
                            }}
                        />

                        <Copyright />
                    </Box>
                </Box>
            </Fade>
        </Collapse>
    )
}