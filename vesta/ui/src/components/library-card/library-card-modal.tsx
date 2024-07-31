'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import ButtonBase from '@mui/material/ButtonBase';
import { useTheme } from '@mui/material/styles';
import Modal from '@mui/material/Modal';
import Typography from '@mui/material/Typography';
import html2canvas from 'html2canvas';

import LibraryCard, { Classes as LibraryCardClasses } from './library-card';
import { TypographyVariant } from '@/lib/utils/types';
import { User } from '../user/models';


export function LibraryCardModal({
    open,
    handleSetOpen,
    user,
}: {
    open: boolean,
    handleSetOpen: (open: boolean) => void,
    user: User,
}) {
    const theme = useTheme()

    const libraryCardRef = React.useRef<HTMLElement>()

    const textTypography: TypographyVariant = 'pBodyMediumWt'

    const handleClickSave = async () => {
        const el = libraryCardRef.current
        if (el) {
            const canvasEl = await html2canvas(el)
            const imageUrl = canvasEl
                .toDataURL('image/png')
                .replace("image/png", "image/octet-stream")

            const filename = `ucsf_data_library_card-${user.name}.png`;

            const a = document.createElement('a');
            a.setAttribute('href', imageUrl);
            a.setAttribute('download', filename);
            a.click();
            a.remove();
        }
    }

    return (
        <Modal
            open={open}
            onClose={(_, reason) => {
                switch (reason) {
                    case 'backdropClick':
                        return
                    case 'escapeKeyDown':
                        handleSetOpen(false)
                        return
                }
            }}
        >
            <Box
                sx={{
                    height: '100%',
                    display: 'flex',
                    flexDirection: 'column',
                    justifyContent: 'center',
                    alignItems: 'center',
                    gap: '11px',
                    [`& .${LibraryCardClasses.root}`]: {
                        position: 'relative',
                        top: open ? '0px' : '100%',
                        transition: theme.transitions.create(
                            'top',
                            {
                                easing: theme.transitions.easing.quint,
                                duration: theme.transitions.duration.quint,
                            },
                        ),
                    },
                }}
            >
                <LibraryCard
                    ref={libraryCardRef}
                    user={user}
                />

                <Box
                    sx={{
                        display: 'flex',
                        flexDirection: 'column',
                        justifyContent: 'center',
                        alignItems: 'center',
                        gap: '10px',
                        px: '10px',
                        py: '10px',
                    }}
                >
                    <ButtonBase
                        onClick={handleClickSave}
                        sx={{
                            width: 'fit-content',
                            px: '10px',
                            py: '8px',
                            borderRadius: '10px',
                            bgcolor: 'ground.grade25',
                        }}
                    >
                        <Typography
                            variant={textTypography}
                            sx={{
                                color: '#FFFDEC',
                            }}
                        >
                            Save as Image
                        </Typography>
                    </ButtonBase>

                    <ButtonBase
                        onClick={() => handleSetOpen(false)}
                        sx={{
                            width: 'fit-content',
                            px: '10px',
                            borderRadius: '10px',
                        }}
                    >
                        <Typography
                            variant={textTypography}
                            sx={{
                                color: '#FFFDEC',
                            }}
                        >
                            Dismiss
                        </Typography>
                    </ButtonBase>
                </Box>
            </Box>
        </Modal>
    )
}