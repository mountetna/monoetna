import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { useTheme } from '@mui/material';
import ButtonBase from '@mui/material/ButtonBase';

import Image, { StaticImageData } from 'next/image';


const iconContainerClass = 'icon-container'
const activationIndicatorClass = 'activation-indicator'


export default function DrawerButton({
    label,
    icon,
    onClick,
    activated,
    open,
}: {
    label: string,
    icon: StaticImageData,
    onClick: () => void,
    activated: boolean,
    open: boolean,
}) {
    const theme = useTheme()

    const transition = theme.transitions.create(
        ['all'],
        {
            easing: theme.transitions.easing.ease,
            duration: theme.transitions.duration.ease,
        }
    )

    // TODO: make this a theme var?
    const height = '52px'

    return (
        <ButtonBase
            className={`drawer-button${activated ? ' activated' : ''}`}
            onClick={onClick}
            sx={{
                display: 'flex',
                flexDirection: 'row',
                minHeight: height,
                maxHeight: height,
                gap: '10px',
                p: '6px 9px',
                pr: '21px',
                bgcolor: 'utilityWhite.main',
                borderRadius: '30px',
                transition,
                '&.activated': {
                    [`& .${iconContainerClass}`]: {
                        bgcolor: 'utilityHighlight.main',
                    },
                    [`& .${activationIndicatorClass}`]: {
                        opacity: 1,
                    },
                },
                '&:hover, &:focus': {
                    bgcolor: 'ground.grade100',
                    [`& .${iconContainerClass}`]: {
                        bgcolor: 'utilityWhiteTransparent25.main',
                    },
                },
            }}
        >
            <Box
                className={iconContainerClass}
                sx={{
                    position: 'relative',
                    display: 'flex',
                    justifyContent: 'center',
                    alignItems: 'center',
                    p: '5px 6.25px',
                    pt: '9px',
                    borderRadius: '50%',
                    transition,
                }}
            >
                <Image
                    src={icon}
                    alt={`Icon for ${label}`}
                    width={20}
                />
                <Box
                    className={activationIndicatorClass}
                    sx={{
                        position: 'absolute',
                        bottom: '1px',
                        right: '1px',
                        opacity: 0,
                        transition,
                    }}
                >
                    <Box
                        sx={{
                            width: '10px',
                            height: '10px',
                            bgcolor: 'utilityWhite.main',
                            borderRadius: '50%',
                        }}
                    />
                    <Box
                        sx={{
                            position: 'absolute',
                            top: '50%',
                            left: '50%',
                            transform: 'translate(-50%, -50%)',
                            width: '6px',
                            height: '6px',
                            bgcolor: 'magenta.grade25',
                            borderRadius: '50%',
                        }}
                    />
                </Box>
            </Box>

            <Typography
                variant='pBody'
                sx={{
                    color: 'ground.grade10',
                }}
            >
                {label}
            </Typography>
        </ButtonBase>
    )
}