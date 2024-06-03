import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';


interface Stat {
    value: string
    label: string
}

export default function SimpleStat({
    primary,
    secondary,
    deltaSign,
    textColor,
}: {
    primary: Stat,
    secondary?: Stat,
    deltaSign?: '+' | '-',
    textColor?: string,
}) {
    return (
        <Box
            className='simple-stat'
            sx={(theme) => ({
                display: 'flex',
                flexDirection: 'column',
                justifyContent: 'space-between',
                p: '16px',
                borderRadius: '30px',
                color: textColor || 'unset',
            })}
        >
            <Typography
                variant='h5'
                sx={(theme) => ({
                    height: '100%'
                    // mb: '1.5em',
                })}
            >
                {primary.label}
            </Typography>
            <Box
                sx={{
                    mb: '5px',
                }}
            >
                <Typography variant='h3Digits'>
                    {primary.value}
                </Typography>
                {deltaSign &&
                    <Box>
                        {deltaSign}
                    </Box>}
                {secondary &&
                    <Box>
                        <span>{secondary.label}</span>
                        <span>{secondary.value}</span>
                    </Box>}
            </Box>
        </Box>
    )
}