import * as React from 'react';
import Box from '@mui/system/Box'
import { Collapse, Fade, Typography, useTheme } from '@mui/material';
import { TransitionProps } from '@mui/material/transitions';
import ButtonBase from '@mui/material/ButtonBase';

import LibraryCardButton, { Classes as LibraryCardButtonClasses } from './library-card-button';
import { User } from '../user/models';
import LibraryCard from './library-card';
import { ClickAwayListener } from '@mui/base';
import { FILE_EXPORT_STATUS, handleExportElementToImage } from '@/lib/utils/file-export';

import downloadSimpleIcon from '/public/images/icons/download-simple.svg'
import Image from 'next/image';


export enum Classes {
    root = 'library-card-tray'
}


interface Props extends React.HTMLAttributes<HTMLDivElement> {
    open: boolean
    onSetOpen: (open: boolean) => void
    user: User | null
    disabled?: boolean,
}


function LibraryCardTray(props: Props, ref: React.ForwardedRef<unknown>) {
    const { open, onSetOpen, user, disabled = false, ...htmlProps } = props
    const theme = useTheme()

    const libraryCardRef = React.useRef<HTMLElement>()
    const [saveImageState, setSaveImageState] = React.useState(FILE_EXPORT_STATUS.idle)

    const handleClickSave = async () => {
        const el = libraryCardRef.current
        if (!el) return

        await handleExportElementToImage(
            el,
            'svg',
            (status) => setSaveImageState(status),
            `UCSF Data Library Card - ${user?.name}`,
        )
    }

    const animationProps: TransitionProps = {
        in: open,
        easing: theme.transitions.easing.quint,
        timeout: theme.transitions.duration.quint,
    }

    const mvmtTransition = theme.transitions.create(
        ['background-color', 'color'],
        {
            easing: theme.transitions.easing.quint,
            duration: theme.transitions.duration.quint,
        },
    )
    const gapTransition = theme.transitions.create(
        ['gap'],
        {
            easing: theme.transitions.easing.ease,
            duration: theme.transitions.duration.ease,
        },
    )
    const opacityTransition = theme.transitions.create(
        ['opacity'],
        {
            easing: theme.transitions.easing.ease,
            duration: theme.transitions.duration.ease,
        },
    )


    return (
        <ClickAwayListener
            onClickAway={() => onSetOpen(false)}
        >
            <Box
                {...htmlProps}
                ref={ref}
                className={Classes.root}
                onClick={() => !open && onSetOpen(true)}
                sx={{
                    display: 'flex',
                    flexDirection: 'column',
                    gap: '13px',
                    pt: '8px',
                    bgcolor: open ? '#EFEFEF' : 'blue.grade50',
                    transition: `${mvmtTransition}, ${gapTransition}`,
                    borderRadius: '0px 0px 20px 20px',
                    cursor: open ? 'unset' : 'pointer',
                    [`.${LibraryCardButtonClasses.root}`]: {
                        bgcolor: open ? '#EFEFEF' : 'blue.grade50',
                        color: open ? 'ground.grade10' : 'utilityHighlight.main',
                        p: '10px',
                        width: '100%',
                        justifyContent: 'flex-start',
                        transition: mvmtTransition,
                    },
                }}
            >
                <Box>
                    {user && (
                        <Collapse
                            {...animationProps}
                        >
                            <Fade
                                {...animationProps}
                            >
                                <Box
                                    sx={{
                                        px: '20px',
                                        py: '10px',
                                    }}
                                >
                                    <LibraryCard
                                        ref={libraryCardRef}
                                        user={user}
                                        variant='2d'
                                    />
                                </Box>
                            </Fade>
                        </Collapse>
                    )}
                </Box>

                <Box
                    sx={{
                        display: 'flex',
                        flexDirection: 'row',
                        justifyContent: 'space-between',
                    }}
                >
                    <LibraryCardButton
                        isLoggedIn={user !== null}
                        onClick={() => onSetOpen(!open)}
                        textOverride={open ? 'Hide your Library Card' : undefined}
                        disabled={disabled}
                    />

                    <ButtonBase
                        onClick={handleClickSave}
                        disabled={!open || saveImageState === FILE_EXPORT_STATUS.inProgress}
                        sx={{
                            px: '16px',
                            py: '8px',
                            gap: '5px',
                            opacity: open ? 1 : 0,
                            transition: opacityTransition,
                        }}
                    >
                        <Image
                            src={downloadSimpleIcon}
                            alt='Download icon'
                            height={22}
                        />

                        <Typography variant='pBody'>
                            Save
                        </Typography>
                    </ButtonBase>
                </Box>
            </Box>
        </ClickAwayListener>
    )
}


export default React.forwardRef(LibraryCardTray)