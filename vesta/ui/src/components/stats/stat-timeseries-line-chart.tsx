import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import * as d3 from 'd3'

import { Instance as StatInstance } from './types';


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
    const dataValues = data.map((datum) => datum.value)
    const dataMin = Math.min(...dataValues)
    const dataMax = Math.max(...dataValues)

    React.useEffect(() => {
        const element = chartContainerRef.current
        if (element === undefined) return

        const width = element.offsetHeight
        const height = element.offsetWidth

        const xScale = d3.scaleTime()
            .domain([data[0].date, data[data.length - 1].date])
            .range([0, width])

        const yScale = d3.scaleLinear()
            .domain([dataMin, dataMax])
            .range([0, height])

        // debugger

    }, [data])

    return (
        <Box
            className='stat-graph'
            sx={(theme) => ({
                display: 'flex',
                flexDirection: 'column',
                p: '16px',
                bgcolor: 'magenta.grade50',
                borderRadius: '30px',
            })}
        >
            <Typography variant='h5'>
                {headingLabel}
            </Typography>
            <Typography variant='h3Digits' component='h3'>
                {headingValue}
            </Typography>
            <Box
                ref={chartContainerRef}
            >

            </Box>
        </Box>
    )
}