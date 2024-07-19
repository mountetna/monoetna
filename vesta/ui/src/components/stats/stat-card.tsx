import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import Image from 'next/image';
import { useTheme } from '@mui/material';

import { useWindowDimensions } from '@/lib/utils/responsive';

import triangleDarkUp from '/public/images/icons/indicator-arrow-dark.svg'
import triangleLightUp from '/public/images/icons/indicator-arrow-light.svg'


interface Stat {
    value: string
    label: string
}

export default function StatCard({
    primary,
    secondary,
    deltaSign,
    deltaColor,
    textColor,
    backgroundColor,
    alwaysShowSecondary = false,
}: {
    primary: Stat,
    secondary?: Stat,
    deltaSign?: '+' | '-',
    deltaColor?: 'light' | 'dark',
    textColor?: string,
    backgroundColor?: string,
    alwaysShowSecondary?: boolean,
}) {
    const theme = useTheme()
    const {
        dimensions: windowDimensions,
        isResizing: isWindowResizing,
    } = useWindowDimensions()

    const transition = theme.transitions.create(
        ['all'],
        {
            duration: theme.transitions.duration.ease,
            easing: theme.transitions.easing.ease,
        },
    )

    const statContainer = React.useRef<HTMLElement>()
    const statPrimaryContainer = React.useRef<HTMLElement>()
    const [
        statPrimaryContainerInitialVerticalOffsetPx,
        setStatPrimaryContainerInitialVerticalOffsetPx
    ] = React.useState<number>(0)

    React.useEffect(() => {
        const _statContainer = statContainer.current
        const _statPrimaryContainer = statPrimaryContainer.current
        if (
            alwaysShowSecondary === true
            || (_statContainer === undefined || _statPrimaryContainer === undefined)
            || [_statContainer.offsetHeight, _statPrimaryContainer.offsetHeight].indexOf(0) >= 0
        ) {
            return
        }

        setStatPrimaryContainerInitialVerticalOffsetPx(
            _statContainer.offsetHeight - _statPrimaryContainer.offsetHeight
        )
    }, [statContainer, statPrimaryContainer, windowDimensions])

    return (
        <Box
            className='stat-card'
            sx={{
                display: 'flex',
                flexDirection: 'column',
                justifyContent: 'space-between',
                p: '16px',
                borderRadius: '30px',
                color: textColor || 'unset',
                bgcolor: backgroundColor || 'unset',
                '& .stat-primary': {
                    position: 'relative',
                    top: `${statPrimaryContainerInitialVerticalOffsetPx}px`,
                    transition: isWindowResizing ? 'unset' : transition,
                },
                '& .delta-icon, & .stat-secondary': {
                    transition: transition,
                    opacity: alwaysShowSecondary ? 1 : 0,
                },
                '&:hover, &:focus': {
                    '& .stat-primary': {
                        position: 'relative',
                        top: 0,
                    },
                    '& .delta-icon, & .stat-secondary': {
                        opacity: 1,
                    },
                },
            }}
        >
            <Typography variant='h5'>
                {primary.label}
            </Typography>
            <Box
                ref={statContainer}
                sx={{
                    mb: '5px',
                }}
            >
                <Box
                    ref={statPrimaryContainer}
                    className='stat stat-primary'
                    sx={{
                        display: 'flex',
                    }}
                >
                    <Typography variant='h3Digits' component='h3'>
                        {primary.value}
                    </Typography>
                    {deltaSign &&
                        <Box
                            className='delta-icon'
                            sx={{
                                display: 'flex',
                                alignItems: 'flex-end',
                                ml: '10px',
                                transform: `rotate(${deltaSign === '+' ? '0' : '180'}deg)`,
                                alignSelf: 'start',
                                width: '24px',
                                height: '24px',
                                '& img': {
                                    width: '100%',
                                    height: 'auto',
                                },
                            }}
                        >
                            <Image
                                src={deltaColor === 'light' ? triangleLightUp : triangleDarkUp}
                                alt={`Triangle icon pointing ${deltaSign === '+' ? 'up' : 'down'}`}
                            />
                        </Box>}
                </Box>
                {secondary &&
                    <Box
                        className='stat stat-secondary'
                    >
                        <Typography variant='h5'>
                            {secondary.value} {secondary.label}
                        </Typography>
                    </Box>}
            </Box>
        </Box>
    )
}