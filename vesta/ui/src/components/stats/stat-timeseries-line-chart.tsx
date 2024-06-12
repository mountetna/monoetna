'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { useTheme } from '@mui/material';
import * as Plot from '@observablehq/plot'

import { Instance as StatInstance } from './types';
import { useWindowDimensions } from '@/lib/utils/responsive';


const dateToLocaleString = (date: Date): string => {
    return [
        date.toLocaleString('local', { month: 'short' }),
        date.toLocaleString('local', { year: 'numeric' }),
    ].join(' ')
}

export default function StatTimeseriesLineChart({
    headingLabel,
    headingValue,
    data,
}: {
    headingLabel: string,
    headingValue: string,
    data: StatInstance<number>[]

}) {
    const chartContainerRef = React.useRef<HTMLElement>()
    const windowDimensions = useWindowDimensions()
    const horizontalPaddingPx = 16

    const theme = useTheme()

    React.useEffect(() => {
        const svgContainer = chartContainerRef.current
        if (svgContainer === undefined) return

        const width = svgContainer.offsetWidth
        const height = svgContainer.offsetHeight

        const plot = Plot.plot({
            width,
            height,
            marginLeft: 0,
            marginRight: 0,
            x: {
                insetLeft: horizontalPaddingPx,
                insetRight: horizontalPaddingPx,
            },
            y: { axis: null },
            marks: [
                Plot.lineY(data, {
                    x: 'date',
                    y: 'value',
                    tip: true,
                    stroke: theme.palette.red.grade25,
                    strokeWidth: 3,
                }),
                Plot.axisX({
                    tickSize: 0,
                    tickFormat: (d: Date) => dateToLocaleString(d),
                }),
                Plot.frame({ anchor: 'bottom' })
            ],
        })

        svgContainer.appendChild(plot)

        return () => plot.remove()
    },
        [data, windowDimensions]
    )

    return (
        <Box
            className='stat-graph'
            sx={(theme) => ({
                display: 'flex',
                flexDirection: 'column',
                // p: `${horizontalPaddingPx}px`,
                bgcolor: 'magenta.grade50',
                borderRadius: '30px',
                overflow: 'hidden',
            })}
        >
            <Box
                sx={{
                    p: `${horizontalPaddingPx}px`,
                }}
            >
                <Typography variant='h5'>
                    {headingLabel}
                </Typography>
                <Typography variant='h3Digits' component='h3'>
                    {headingValue}
                </Typography>
            </Box>
            <Box
                ref={chartContainerRef}
                sx={{
                    flexGrow: 1,
                    overflow: 'hidden',
                }}
            />
        </Box>
    )
}