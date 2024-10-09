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


export default function DrawerNav({
    accessUrl,
    open,
    onClose,
    sx = {},
    onClickNavLink,
    onSetLibraryCardModalOpen,
}: {
    accessUrl: string,
    open: boolean,
    onClose: () => void,
    sx?: SxProps,
    onClickNavLink?: () => void,
    onSetLibraryCardModalOpen?: (open: boolean) => void,
}) {
    const theme = useTheme()
    const user = useUser()

    const [libraryCardModalOpen, setLibraryCardModalOpen] = React.useState(false)
    const handleSetLibraryCardModalOpen = (open: boolean) => {
        setLibraryCardModalOpen(open)
        onSetLibraryCardModalOpen && onSetLibraryCardModalOpen(open)
    }

    const animationProps: TransitionProps = {
        in: open,
        easing: theme.transitions.easing.quint,
        timeout: theme.transitions.duration.quint,
    }

    return (
        <Collapse
            component='nav'
            aria-label='Mobile'
            // open={open}
            // onClose={onClose}
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
                        [theme.breakpoints.up('tablet')]: {
                            pb: '16px',
                        },
                        ...sx,
                    }}
                >
                    <Box
                        sx={{
                            display: 'flex',
                            flexDirection: 'column',
                            gap: '24px',
                            px: '16px',
                            [`& .${LibraryCardButtonClasses.root}`]: {
                                width: 'fit-content',
                            },
                        }}
                    >
                        <DLNav
                            linkTypography='h2'
                            onClickNavLink={onClickNavLink}
                            sx={{
                                display: 'flex',
                                flexDirection: 'column',
                                alignItems: 'start',
                                gap: '16px',
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
                            }}
                        />

                        <Box
                            sx={{
                                [theme.breakpoints.up('tablet')]: {
                                    display: 'none',
                                },
                            }}
                        >
                            <LibraryCardButton
                                accessUrl={accessUrl}
                                isLoggedIn={user !== null}
                                onClick={() => handleSetLibraryCardModalOpen(!libraryCardModalOpen)}
                            />

                            {user && (
                                <LibraryCardModal
                                    open={libraryCardModalOpen}
                                    handleSetOpen={handleSetLibraryCardModalOpen}
                                    user={user}
                                />
                            )}
                        </Box>
                    </Box>

                    <Box
                        sx={{
                            display: 'flex',
                            flexDirection: 'column',
                            gap: '16px',
                            px: '16px',
                            [theme.breakpoints.up('tablet')]: {
                                px: '0px'
                            },
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