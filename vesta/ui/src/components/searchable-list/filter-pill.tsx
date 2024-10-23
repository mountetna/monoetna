'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { useTheme } from '@mui/material';
import ButtonBase from '@mui/material/ButtonBase';
import Image from 'next/image';

import xDark from '/public/images/icons/x-dark.svg'


const removeButtonClass = 'remove-button'


export default function FilterPill({
    label,
    removeable = true,
    onClickRemove,
}: {
    label: string,
    removeable?: boolean,
    onClickRemove?: () => void,
}) {
    const theme = useTheme()

    const transition = theme.transitions.create(
        ['all'],
        {
            easing: theme.transitions.easing.ease,
            duration: theme.transitions.duration.ease,
        }
    )

    return (
        <ButtonBase
            onClick={() => onClickRemove && onClickRemove()}
            sx={{
                display: 'flex',
                flexDirection: 'row',
                gap: '12px',
                px: '14px',
                py: '12px',
                pr: '21px',
                bgcolor: 'utilityWhite.main',
                borderRadius: '30px',
                transition,
                '&:hover, &:focus': {
                    bgcolor: 'ground.grade100',
                    [`& .${removeButtonClass}`]: {
                        bgcolor: 'utilityWhiteTransparent25.main',
                    },
                },
            }}
        >
            {removeable && <Box
                className={removeButtonClass}
                sx={{
                    position: 'relative',
                    display: 'flex',
                    justifyContent: 'center',
                    alignItems: 'center',
                    p: '8px',
                    borderRadius: '50%',
                    transition,
                }}
            >
                <Image
                    src={xDark}
                    alt={`"X" icon`}
                    width={8}
                    height={8}
                />
            </Box>}

            <Typography
                variant='pBody'
                sx={{
                    color: 'ground.grade10',
                }}
            >
                {label}
            </Typography>
        </ButtonBase >
    )
}