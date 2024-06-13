'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { alpha, useMediaQuery, useTheme } from '@mui/material';
import { useParentSize } from '@visx/responsive'
import {
    Axis,
    AnimatedLineSeries,
    XYChart,
    Tooltip,
} from '@visx/xychart';

import { Instance as StatInstance } from './types';
import { DEFAULT_LETTER_SPACING } from '@/theme';
import Dropdown from '../inputs/dropdown';


const dateToLocaleString = (date: Date, includeDay: boolean = true): string => {
    const parts = [
        date.toLocaleString('local', { month: 'short' }),
        date.toLocaleString('local', { year: 'numeric' }),
    ]

    if (includeDay) {
        parts.unshift(
            date.toLocaleString('local', { day: 'numeric' })
        )
    }

    return parts.join(' ')
}

const dayInMs = 24 * 60 * 60 * 1000


export default function StatTimeseriesLineChart({
    headingLabel,
    headingValue,
    data,
}: {
    headingLabel: string,
    headingValue: string,
    data: StatInstance<number>[]

}) {
    const [dropwdownVal, setDropdownVal] = React.useState<keyof typeof timeWindowOptions>('all')

    const {
        parentRef: chartContainerRef,
        width: chartContainerWidth,
        height: chartContainerHeight,
    } = useParentSize({ debounceTime: 100, })
    const paddingPx = 16

    const theme = useTheme()
    const isMobile = useMediaQuery(theme.breakpoints.down('tablet'))
    const isTablet = useMediaQuery(theme.breakpoints.between(
        theme.breakpoints.values.tablet,
        theme.breakpoints.values.desktop,
    ))
    const isDesktop = useMediaQuery(theme.breakpoints.between(
        theme.breakpoints.values.desktop,
        theme.breakpoints.values.desktopLg,
    ))
    // const isDesktopLg = useMediaQuery(theme.breakpoints.up('desktopLg'))

    const accessors = {
        xAccessor: (d: StatInstance<number>) => d.date,
        yAccessor: (d: StatInstance<number>) => d.value,
    };

    const timeWindowOptions = {
        week: { short: 'Week', full: 'Week', timeWindowStart: new Date(Date.now() - 7 * dayInMs) },
        month: { short: 'Month', full: 'Month', timeWindowStart: new Date(Date.now() - 30 * dayInMs) },
        year: { short: 'Year', full: 'Year', timeWindowStart: new Date(Date.now() - 365 * dayInMs) },
        all: { short: 'All', full: 'All time', timeWindowStart: new Date(0) },
    }

    const filteredData = data.filter((val) => (
        val.date >= timeWindowOptions[dropwdownVal].timeWindowStart
    ))

    const showDayInAxis = ['month', 'week'].indexOf(dropwdownVal) >= 0 || filteredData.length <= 30

    return (
        <Box
            className='stat-graph'
            sx={(theme) => ({
                display: 'flex',
                flexDirection: 'column',
                bgcolor: 'magenta.grade50',
                borderRadius: '30px',
                overflow: 'hidden',
                [theme.breakpoints.up('desktop')]: {
                    pb: `${paddingPx}px`,
                },
            })}
        >
            <Box
                sx={{
                    mb: '10px',
                    p: `${paddingPx}px`,
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
            >
                <XYChart
                    width={chartContainerWidth}
                    height={chartContainerHeight}
                    xScale={{ type: 'time' }}
                    yScale={{ type: 'linear' }}
                    margin={{ top: 0, right: 0, bottom: 32, left: 0 }}
                >
                    <Axis
                        orientation='bottom'
                        hideZero
                        hideTicks
                        stroke={alpha('#000', 0.15)}
                        tickFormat={(v: Date) => dateToLocaleString(v, showDayInAxis)}
                        numTicks={isMobile ? 3 : isTablet ? 3 : isDesktop ? 6 : 7}
                        tickLabelProps={{
                            style: {
                                fontFamily: theme.typography.pLarge.fontFamily,
                                fontSize: isMobile ? 18 : isTablet ? 20 : 22,
                                letterSpacing: DEFAULT_LETTER_SPACING,
                                lineHeight: isMobile ? '150%' : isTablet ? '148%' : '148%',
                            },
                            stroke: theme.palette.ground.grade10,
                            textAnchor: 'start',
                            dy: 16,
                        }}

                    />
                    <AnimatedLineSeries
                        dataKey='Users'
                        data={filteredData}
                        stroke={theme.palette.red.grade25}
                        strokeWidth={3}
                        {...accessors}
                    />
                    <Tooltip<StatInstance<number>>
                        renderTooltip={({ tooltipData }) => {
                            if (tooltipData?.nearestDatum?.datum === undefined) return <div></div>

                            return (
                                <div>
                                    {accessors.yAccessor(tooltipData.nearestDatum.datum)} users
                                    {', '}
                                    {accessors.xAccessor(tooltipData.nearestDatum.datum).toLocaleDateString()}
                                </div>
                            )
                        }}
                    />
                </XYChart>
            </Box>
            <Box
                sx={(theme) => ({
                    p: `${paddingPx}px`,
                    [theme.breakpoints.up('desktop')]: {
                        display: 'none',
                    },
                    '& .dropdown-button-content': {
                        justifyContent: 'center',
                    },
                })}
            >
                <Dropdown
                    options={
                        Object.entries(timeWindowOptions).map(
                            ([k, v]) => ({ value: k, label: v.full })
                        )
                    }
                    value={dropwdownVal}
                    onChange={(newVal) => setDropdownVal(newVal)}
                />
            </Box>
        </Box>
    )
}