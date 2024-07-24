import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { useTheme } from '@mui/material';
import ButtonBase from '@mui/material/ButtonBase';

import Image, { StaticImageData } from 'next/image';


export enum Class {
    buttonBase = 'drawer-button-base',
    buttonActivated = 'drawer-button-activated',
    buttonOpen = 'drawer-button-apen',
    iconContainer = 'icon-container',
    activationIndicator = 'activation-indicator',
    activationIndicatorBorder = 'activation-indicator-border',
}


export default function DrawerButton({
    label,
    iconLight,
    iconDark,
    onClick,
    activated,
    open,
}: {
    label: string,
    iconLight: StaticImageData,
    iconDark: StaticImageData,
    onClick: () => void,
    activated: boolean,
    open: boolean,
}) {
    const theme = useTheme()

    const transition = theme.transitions.create(
        'all',
        {
            easing: theme.transitions.easing.ease,
            duration: theme.transitions.duration.ease,
        }
    )

    const classes = [Class.buttonBase]
    if (activated) {
        classes.push(Class.buttonActivated)
    }
    if (open) {
        classes.push(Class.buttonOpen)
    }

    // TODO: make this a theme var?
    const height = '52px'

    const iconAlt = `${open ? 'Light' : 'Dark'} icon for ${label}`

    return (
        <ButtonBase
            className={classes.join(' ')}
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
                [`&.${Class.buttonActivated}`]: {
                    [`& .${Class.iconContainer}`]: {
                        bgcolor: 'utilityHighlight.main',
                    },
                    [`& .${Class.activationIndicator}`]: {
                        opacity: 1,
                    },
                },
                '&:hover, &:focus': {
                    [`&:not(.${Class.buttonOpen})`]: {
                        bgcolor: 'ground.grade100',
                        [`& .${Class.iconContainer}`]: {
                            bgcolor: 'utilityWhiteTransparent25.main',
                        },
                    },
                },
                [`&.${Class.buttonOpen}`]: {
                    bgcolor: 'ground.grade10',
                    [`& .${Class.iconContainer}`]: {
                        bgcolor: 'ground.grade25',
                    },
                    [`& .${Class.activationIndicatorBorder}`]: {
                        bgcolor: 'ground.grade10',
                    },
                    '& *': {
                        color: 'utilityHighlight.main',
                    },
                },
            }}
        >
            <Box
                className={Class.iconContainer}
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
                <Box
                    sx={{
                        position: 'relative',
                    }}
                >
                    <Image
                        src={iconDark}
                        alt={iconAlt}
                        width={20}
                    />

                    <Box
                        sx={{
                            position: 'absolute',
                            top: 0,
                            left: 0,
                            opacity: open ? 1 : 0,
                            transition,
                        }}
                    >
                        <Image
                            src={iconLight}
                            alt={iconAlt}
                            width={20}
                        />
                    </Box>
                </Box>
                <Box
                    className={Class.activationIndicator}
                    sx={{
                        position: 'absolute',
                        bottom: '1px',
                        right: '1px',
                        opacity: 0,
                        transition,
                    }}
                >
                    <Box
                        className={Class.activationIndicatorBorder}
                        sx={{
                            width: '10px',
                            height: '10px',
                            bgcolor: 'utilityWhite.main',
                            borderRadius: '50%',
                            transition,
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