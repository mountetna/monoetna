'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { Collapse, Fade, useTheme } from '@mui/material';


export default function ProjectStatus({
    currentStatus,
    allStatuses,
    variant,
    visible,
}: {
    currentStatus: string,
    allStatuses: string[],
    variant: 'default' | 'compact'
    visible: boolean,
}) {
    const theme = useTheme()

    const currentStatusIdx = allStatuses.indexOf(currentStatus)

    let mainContent
    switch (variant) {
        case 'default':
            mainContent = (
                <Box
                    sx={{
                        display: 'grid',
                        gridTemplateColumns: 'repeat(12, 1fr)',
                        '& .stage': {
                            display: 'flex',
                            flexDirection: 'column',
                            gap: '13.35px',
                            gridColumn: 'span 3',
                        },
                        '& .bubble-container': {
                            position: 'relative',
                            width: '100%',
                            height: '19px',
                        },
                        '& .bubble': {
                            width: '100%',
                            height: '100%',
                            bgcolor: 'ground.grade50',
                        },
                        '& .bubble-overlay': {
                            position: 'absolute',
                            left: 0,
                            top: 0,
                            height: '100%',
                            bgcolor: 'ground.grade10',
                        },
                        '& .bubble, & .bubble-overlay': {
                            borderRadius: '40px',
                        },
                        '& .stage.passed .bubble-overlay, & .stage.current:last-child .bubble-overlay': {
                            width: '100%',
                        },
                        '& .stage.current:first-child .bubble-overlay': {
                            width: '10%',
                        },
                        '& .tick-text-container': {
                            display: 'flex',
                            flexDirection: 'column',
                            gap: '7.19px',
                        },
                        '& .tick': {
                            width: '1px',
                            height: '12px',
                            bgcolor: 'ground.grade25',
                            opacity: 0.4,
                        },
                        '& .text': {
                            color: 'ground.grade25',
                            opacity: 0.4,
                        },
                        '& .stage.current': {
                            '& .tick': {
                                bgcolor: 'ground.grade10',
                                opacity: 1,
                            },
                            '& .text': {
                                color: 'ground.grade10',
                                opacity: 1,
                            },
                        },
                    }}
                >
                    {allStatuses.map((currentStatus, i) => (
                        <Box
                            key={i}
                            className={`stage${currentStatusIdx > i ? ' passed' : currentStatusIdx === i ? ' current' : ''}`}
                        >
                            <Box className='bubble-container'>
                                <Box className='bubble'></Box>
                                <Fade
                                    in={visible}
                                    easing={theme.transitions.easing.ease}
                                    timeout={theme.transitions.duration.ease}
                                    style={{
                                        transitionDelay: `${theme.transitions.duration.ease / allStatuses.length * i + theme.transitions.duration.ease / 1.1}ms`,
                                    }}
                                >
                                    <Box className='bubble-overlay'></Box>
                                </Fade>
                            </Box>

                            <Box className='tick-text-container'>
                                <Box className='tick'></Box>

                                <Typography
                                    variant='pMediumBoldWt'
                                    className='text'
                                >
                                    {currentStatus}
                                </Typography>
                            </Box>
                        </Box>
                    ))}
                </Box>
            )
            break

        case 'compact':
            mainContent = (
                <Box
                    sx={{
                        display: 'flex',
                        flexDirection: 'row',
                        gap: '10px',
                    }}
                >
                    <Typography
                        variant='pLarge'
                        sx={{
                            color: 'ground.grade100',
                            bgcolor: 'ground.grade25',
                            borderRadius: '30px',
                            px: '9px',
                        }}
                    >
                        {`${currentStatusIdx + 1} of ${allStatuses.length}`}
                    </Typography>

                    <Typography
                        variant='pMediumMediumWt'
                    >
                        {currentStatus}
                    </Typography>
                </Box>
            )
            break
    }

    return (
        <Box
            sx={{
                display: 'flex',
                flexDirection: 'column',
                gap: '8px',
                pb: '24px',
                [theme.breakpoints.up('tablet')]: {
                    gap: '16.43px',
                },
                [theme.breakpoints.up('desktop')]: {
                    pb: '75px',
                },
            }}
        >
            <Typography
                variant='pMediumBoldWt'
                component='div'
            >
                Project Status
            </Typography>

            {mainContent}
        </Box>
    )
}