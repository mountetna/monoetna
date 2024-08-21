'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { alpha, Breakpoint, Fade, useMediaQuery, useTheme } from '@mui/material';
import ButtonBase from '@mui/material/ButtonBase';
import { useParentSize } from '@visx/responsive'
import { AnimatedLineSeries, XYChart, Tooltip, TooltipProvider } from '@visx/xychart';
import { scaleTime } from '@visx/scale'
import { ScaleTime } from '@visx/vendor/d3-scale'
import Image from 'next/image';

import { Instance as StatInstance } from './types';
import Dropdown from '../inputs/dropdown';
import indicatorArrowDark from '/public/images/icons/indicator-arrow-dark.svg'
import TooltipContent from '../tooltip/tooltip-content';
import { useBreakpoint } from '@/lib/utils/responsive';


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

const IndicatorArrowDark = (
    <Box
        sx={{
            display: 'flex',
            transform: 'rotate(180deg)',
            '& img': {
                width: '19px',
                height: 'auto',
            },
        }}
        component='span'
    >
        <Image
            src={indicatorArrowDark}
            alt='Indicator arrow pointing down'
        />
    </Box>
)


function TimeAxis({
    scale,
    numTicks,
    showDay,
}: {
    scale: ScaleTime<number, number, never>,
    numTicks: number,
    showDay: boolean,
}) {
    return (
        <Box
            sx={{
                display: 'flex',
                justifyContent: 'space-between',
                textAlign: 'center',
            }}
        >
            {Array(numTicks).fill(null).map((_, index) => {
                const date = scale.invert(index / (numTicks - 1) * scale.range()[1])
                const asStr = dateToLocaleString(date, showDay)

                return (
                    <Typography
                        key={index}
                        variant='pLarge'
                        sx={{
                            mx: '8px',
                            '&:first-child': {
                                ml: '0',
                            },
                            '&:last-child': {
                                mr: '0',
                            },
                        }}
                    >
                        {asStr}
                    </Typography>
                )
            })}
        </Box>
    )
}


interface StatTimeseriesLineChartProps {
    data: StatInstance<number>[];
    dataLabelSingular: string;
    dataLabelPlural: string;
    headingLabel: string;
    headingValue: number;
}


function _StatTimeseriesLineChart({
    data,
    dataLabelSingular,
    dataLabelPlural,
    headingLabel,
    headingValue,
}: StatTimeseriesLineChartProps) {
    const [timeWindowVal, setTimeWindowVal] = React.useState<keyof typeof timeWindowOptions>('all')

    const theme = useTheme()
    const transition = theme.transitions.create(
        ['all'],
        {
            duration: theme.transitions.duration.ease,
            easing: theme.transitions.easing.ease,
        },
    )

    const {
        parentRef: chartContainerRef,
        width: chartContainerWidth,
        height: chartContainerHeight,
    } = useParentSize({ debounceTime: 500, })
    const paddingPx = 16

    const isMobile = useMediaQuery(theme.breakpoints.down('tablet'))
    const isTablet = useMediaQuery(theme.breakpoints.between(
        theme.breakpoints.values.tablet,
        theme.breakpoints.values.desktop,
    ))
    const isDesktop = useMediaQuery(theme.breakpoints.between(
        theme.breakpoints.values.desktop,
        theme.breakpoints.values.desktopLg,
    ))

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
        val.date >= timeWindowOptions[timeWindowVal].timeWindowStart
    ))

    const xScale = scaleTime({
        domain: [filteredData[0].date, filteredData[filteredData.length - 1].date],
        range: [0, chartContainerWidth],
    })

    const showDayInAxis = ['month', 'week'].indexOf(timeWindowVal) >= 0 || filteredData.length <= 30

    const numTicksWithDay = isMobile ? 3 : isTablet ? 3 : isDesktop ? 6 : 7
    const numTicksWODay = isMobile ? 3 : isTablet ? 3 : isDesktop ? 8 : 8

    return (
        <Box
            className='stat-graph'
            sx={{
                display: 'flex',
                flexDirection: 'column',
                bgcolor: 'magenta.grade50',
                borderRadius: '30px',
                overflow: 'hidden',
            }}
        >
            <Box
                sx={{
                    display: 'flex',
                    justifyContent: 'space-between',
                    alignItems: 'flex-start',
                    mb: '10px',
                    p: `${paddingPx}px`,
                }}
            >
                <Box
                    sx={{
                        color: 'ground.grade10',
                    }}
                >
                    <Typography variant='h5'>
                        {headingLabel}
                    </Typography>
                    <Typography variant='h3Digits' component='h3'>
                        {headingValue} {headingValue === 1 ? dataLabelSingular : dataLabelPlural}
                    </Typography>
                </Box>
                <Box
                    sx={(theme) => ({
                        display: 'none',
                        [theme.breakpoints.up('desktop')]: {
                            display: 'block',
                        },
                    })}
                    role='listbox'
                >
                    {Object.entries(timeWindowOptions).map(([k, v]) => {
                        const isSelected = k === timeWindowVal

                        return (
                            // @ts-ignore
                            <ButtonBase
                                key={k}
                                className={`time-window-option time-window-option-desktop${isSelected ? ' selected' : ''}`}
                                sx={{
                                    color: 'ground.grade10',
                                    px: '20px',
                                    '&:not(:last-child)': {
                                        mr: '10px',
                                    },
                                    borderRadius: '40px',
                                    transition: transition,
                                    '&.selected': {
                                        bgcolor: 'utilityWhite.main',
                                    },
                                }}
                                onClick={_ => {
                                    // @ts-ignore
                                    setTimeWindowVal(k)
                                }}
                                role='option'
                                aria-selected={isSelected}
                            >
                                <Typography variant='pLarge'>
                                    {v.short}
                                </Typography>
                            </ButtonBase>
                        )
                    })}
                </Box>
            </Box>
            <Box
                ref={chartContainerRef}
                sx={{
                    flexGrow: 1,
                    overflow: 'hidden',
                    '& .visx-tooltip': {
                        position: 'absolute',
                    },
                }}
            >
                {chartContainerWidth > 0 && chartContainerHeight > 0 &&
                    <XYChart
                        width={chartContainerWidth}
                        height={chartContainerHeight}
                        xScale={{ type: 'time' }}
                        yScale={{ type: 'linear' }}
                        margin={{ top: 0, right: 0, bottom: -45, left: 0 }}
                    >
                        <AnimatedLineSeries
                            dataKey='Users'
                            data={filteredData}
                            stroke={theme.palette.red.grade25}
                            strokeWidth={3}
                            {...accessors}
                        />
                        <Tooltip<StatInstance<number>>
                            unstyled
                            applyPositionStyle
                            offsetTop={-55}
                            showDatumGlyph
                            renderTooltip={({ tooltipData }) => {
                                if (tooltipData?.nearestDatum?.datum === undefined) return <div></div>
                                const xVal = accessors.xAccessor(tooltipData.nearestDatum.datum)
                                const yVal = accessors.yAccessor(tooltipData.nearestDatum.datum)

                                // TODO: figure out how to fade out
                                return (
                                    <Fade
                                        in={true}
                                        easing={theme.transitions.easing.ease}
                                        timeout={theme.transitions.duration.ease}
                                    >
                                        <Box>
                                            <TooltipContent
                                                variant='simple'
                                            >
                                                <Box>
                                                    {yVal} {yVal === 1 ? dataLabelSingular : dataLabelPlural}
                                                    {', '}
                                                    {xVal.toLocaleDateString()}
                                                </Box>
                                            </TooltipContent>
                                        </Box>
                                    </Fade>
                                )
                            }}

                        />
                    </XYChart>}
            </Box>
            {/* TODO: replace with visx or d3 axis? */}
            <Box
                sx={{
                    borderTop: `1px solid ${alpha('#000', 0.15)}`,
                    mt: '10px',
                    p: `${paddingPx}px`,
                }}
            >
                {/* Prevent rendering when xScale's range is [0, 0] */}
                {chartContainerWidth > 0 && <TimeAxis
                    scale={xScale}
                    numTicks={showDayInAxis ? numTicksWithDay : numTicksWODay}
                    showDay={showDayInAxis}
                />}
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
                    value={timeWindowVal}
                    onChange={(newVal) => setTimeWindowVal(newVal)}
                    icon={IndicatorArrowDark}
                    listboxStyles={{
                        backgroundColor: theme.palette.utilityHighlight.main,
                        border: 'none',
                    }}
                    optionStyles={{
                        textAlign: 'center',
                    }}
                    optionTypographyVariant='pSubtitleBoldWt'
                    selectedOptionTypographyVariant='pSubtitle'
                />
            </Box>
        </Box>
    )
}

export default function StatTimeseriesLineChart(props: StatTimeseriesLineChartProps) {
    const breakpoint = useBreakpoint()
    const isDesktop = (['desktop', 'desktopLg'] as Breakpoint[]).includes(breakpoint)

    return (
        <TooltipProvider hideTooltipDebounceMs={isDesktop ? 0 : 2000}>
            <_StatTimeseriesLineChart {...props} />
        </TooltipProvider>
    )
}