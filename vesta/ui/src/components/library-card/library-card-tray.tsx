import * as React from 'react';
import Box from '@mui/system/Box'
import { Collapse, Fade, useTheme } from '@mui/material';
import { TransitionProps } from '@mui/material/transitions';

import LibraryCardButton, { Classes as LibraryCardButtonClasses } from './library-card-button';
import { User } from '../user/models';
import LibraryCard from './library-card';
import { ClickAwayListener } from '@mui/base';


export enum Classes {
    root = 'library-card-tray'
}


export default function LibraryCardTray({
    user,
}: {
    user: User | null
}) {
    const [open, setOpen] = React.useState(false)

    const theme = useTheme()
    const animationProps: TransitionProps = {
        in: open,
        easing: theme.transitions.easing.quint,
        timeout: theme.transitions.duration.quint,
    }

    const transition = theme.transitions.create(
        ['background-color', 'color'],
        {
            easing: theme.transitions.easing.quint,
            duration: theme.transitions.duration.quint,
        },
    )

    return (
        <ClickAwayListener
            onClickAway={() => setOpen(false)}
        >
            <Box
                className={Classes.root}
                sx={{
                    position: 'relative',
                    display: 'flex',
                    flexDirection: 'column',
                    bgcolor: open ? '#EFEFEF' : 'blue.grade50',
                    transition,
                    borderRadius: '0px 0px 20px 20px',
                    [`.${LibraryCardButtonClasses.root}`]: {
                        bgcolor: open ? '#EFEFEF' : 'blue.grade50',
                        color: open ? 'ground.grade10' : 'utilityHighlight.main',
                        transition,
                    },
                }}
            >
                <Box
                    sx={{
                        // overflow: 'visible',
                        p: '8px',
                    }}
                >
                    {user && (
                        <Collapse
                            {...animationProps}
                        >
                            <Fade
                                {...animationProps}
                            >
                                <Box>

                                    <LibraryCard
                                        user={user}
                                    />
                                </Box>
                            </Fade>
                        </Collapse>
                    )}
                </Box>

                <LibraryCardButton
                    isLoggedIn={user !== null}
                    onClick={() => setOpen(!open)}
                    textOverride={open ? 'Close' : undefined}
                />
            </Box>
        </ClickAwayListener>
    )
}