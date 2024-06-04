import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import Image from 'next/image';

import carrotDarkUp from '/public/images/icons/indicator-arrow-dark.svg'
import carrotLightUp from '/public/images/icons/indicator-arrow-light.svg'


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
}: {
    primary: Stat,
    secondary?: Stat,
    deltaSign?: '+' | '-',
    deltaColor?: 'light' | 'dark',
    textColor?: string,
    backgroundColor?: string,
}) {
    return (
        <Box
            className='simple-stat'
            sx={{
                display: 'flex',
                flexDirection: 'column',
                justifyContent: 'space-between',
                p: '16px',
                borderRadius: '30px',
                color: textColor || 'unset',
                bgcolor: backgroundColor || 'unset',
            }}
        >
            <Typography variant='h5'>
                {primary.label}
            </Typography>
            <Box
                sx={{
                    mb: '5px',
                }}
            >
                <Box
                    sx={{
                        display: 'flex',
                    }}
                >
                    <Typography variant='h3Digits'>
                        {primary.value}
                    </Typography>
                    {deltaSign &&
                        <Box
                            sx={{
                                ml: '10px',
                                transform: `rotate(${deltaSign === '+' ? '0' : '180'}deg)`
                            }}
                        >
                            <Image
                                src={deltaColor === 'light' ? carrotLightUp : carrotDarkUp}
                                alt={`Carrot icon pointing ${deltaSign === '+' ? 'up' : 'down'}`}
                            />
                        </Box>}
                </Box>
                {secondary &&
                    <Box>
                        <Typography variant='h5'>
                            {secondary.value} {secondary.label}
                        </Typography>
                    </Box>}
            </Box>
        </Box>
    )
}