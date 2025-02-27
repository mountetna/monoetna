'use client'

import * as React from 'react'
import { styled } from '@mui/material';
import Box from '@mui/system/Box'
import ButtonBase from '@mui/material/ButtonBase';
import Fade from '@mui/material/Fade';
import Slide from '@mui/material/Slide';
import { useTheme } from '@mui/material/styles';
import { Modal as BaseModal } from '@mui/base/Modal';
import Typography from '@mui/material/Typography';

import LibraryCard from './library-card';
import { TypographyVariant } from '@/lib/utils/types';
import { User } from '../user/models';
import { FILE_EXPORT_STATUS, handleExportElementToImage } from '@/lib/utils/file-export';


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

    const textTypography: TypographyVariant = 'pBodyMediumWt'
    
    const libraryCardRef = React.useRef<HTMLElement>()
    const [saveImageState, setSaveImageState] = React.useState(FILE_EXPORT_STATUS.idle)

    const handleClickSave = async () => {
        const el = libraryCardRef.current
        if (!el) return

        await handleExportElementToImage(
            el,
            'svg',
            (status) => setSaveImageState(status),
            `UCSF Data Library Card - ${user.name}`,
        )
    }

    return (
        <Modal
            aria-labelledby="Data Library Card Display"
            aria-describedby="Shows your Data Library Card and with the ability to save an image"
            open={open}
            onClose={() => handleSetOpen(false)}
            closeAfterTransition
            slots={{ backdrop: StyledBackdrop }}
        >
            <Fade
                in={open}
                easing={theme.transitions.easing.ease}
                timeout={theme.transitions.duration.quint}
            >
                <Box
                    sx={{
                        height: '100%',
                        display: 'flex',
                        flexDirection: 'column',
                        justifyContent: 'center',
                        alignItems: 'center',
                        gap: '11px',
                    }}
                >
                    <Slide
                        in={open}
                        direction='up'
                        easing={theme.transitions.easing.quint}
                        timeout={theme.transitions.duration.quint}
                    >
                        <Box>
                            <LibraryCard
                                ref={libraryCardRef}
                                user={user}
                                variant='3d'
                            />
                        </Box>
                    </Slide>

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
            </Fade>
        </Modal>
    )
}

const Backdrop = React.forwardRef<HTMLDivElement, { open?: boolean }>(
    (props, ref) => {
        const { open, ...other } = props;

        const theme = useTheme()

        return (
            <Fade
                in={open}
                easing={theme.transitions.easing.ease}
                timeout={theme.transitions.duration.quint}
            >
                <div ref={ref} {...other} />
            </Fade>
        );
    },
);
Backdrop.displayName = 'Backdrop';

const Modal = styled(BaseModal)`
  position: fixed;
  z-index: 1300;
  inset: 0;
  display: flex;
  align-items: center;
  justify-content: center;
`;

const StyledBackdrop = styled(Backdrop)`
  z-index: -1;
  position: fixed;
  inset: 0;
  background-color: rgba(0, 0, 0, 0.7);
  -webkit-tap-highlight-color: transparent;
`;