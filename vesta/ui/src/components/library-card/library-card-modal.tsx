'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import ButtonBase from '@mui/material/ButtonBase';
import { useTheme } from '@mui/material/styles';
import Modal from '@mui/material/Modal';
import Typography from '@mui/material/Typography';
import { toPng } from 'html-to-image';

import LibraryCard, { Classes as LibraryCardClasses } from './library-card';
import { TypographyVariant } from '@/lib/utils/types';
import { User } from '../user/models';
import { FILE_EXPORT_STATUS } from '@/lib/utils/file-export';


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
    const [saveImageState, setSaveImageState] = React.useState(FILE_EXPORT_STATUS.idle)

    const textTypography: TypographyVariant = 'pBodyMediumWt'

    const handleClickSave = async () => {
        setSaveImageState(FILE_EXPORT_STATUS.inProgress)

        const el = libraryCardRef.current
        if (!el) return

        try {
            const imageUrl = await toPng(el)
            const filename = `ucsf_data_library_card-${user.name}.png`;

            const a = document.createElement('a');
            a.setAttribute('href', imageUrl);
            a.setAttribute('download', filename);
            a.click();
            a.remove();

            setSaveImageState(FILE_EXPORT_STATUS.success)
        } catch (error) {
            setSaveImageState(FILE_EXPORT_STATUS.error)
            console.error('Error saving image', error)
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
                        disabled={saveImageState === FILE_EXPORT_STATUS.inProgress}
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