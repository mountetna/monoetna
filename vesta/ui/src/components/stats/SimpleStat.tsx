import * as React from 'react'
import Box from '@mui/system/Box'


interface Stat {
    value: string
    label: string
}

export default function SimpleStat({
    primary,
    secondary,
    deltaSign,
}: {
    primary: Stat,
    secondary?: Stat,
    deltaSign?: '+' | '-',
}) {
    return (
        <Box>
            <Box>
                <span>{primary.label}</span>
                <span>{primary.value}</span>
            </Box>
            {deltaSign && <Box>
                {deltaSign}
            </Box>}
            {secondary &&
                <Box>
                    <span>{secondary.label}</span>
                    <span>{secondary.value}</span>
                </Box>}
        </Box>
    )
}