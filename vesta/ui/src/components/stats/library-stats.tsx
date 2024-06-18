'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Container from '@mui/system/Container'
import Typography from '@mui/material/Typography';

import StatCard from './stat-card'
import StatsCarousel, { Stats as StatsProp } from './stats-carousel'
import StatTimeseriesLineChart from './stat-timeseries-line-chart'
import { Instance } from './types';
import theme, { headerMargins } from '@/theme';
import { SIValue, roundValueToNearestSIPrefix } from '@/lib/utils/units';
import ThemeProjectBreakdownChart, { ThemeData } from './theme-project-breakdown-chart';


export type StatsTimeseries = Record<keyof StatsProp, Instance<number>[]>

export default function LibraryStats({
    stats,
    themeProjectBreakdown,
}: {
    stats: StatsTimeseries,
    themeProjectBreakdown: ThemeData[],
}) {
    const darkText = theme.palette.ground.grade10
    const lightText = theme.palette.utilityHighlight.main
    const itemGap = '16px'

    const latestStats = {} as Record<keyof StatsProp, SIValue>
    const since7DaysAgo = {} as Record<keyof StatsProp, SIValue>
    const carouselStats = {} as Record<keyof StatsProp, number>

    for (const [k, v] of (Object.entries(stats) as [keyof StatsProp, Instance<number>[]][])) {

        const latest = v[v.length - 1].value
        latestStats[k] = roundValueToNearestSIPrefix(latest)
        since7DaysAgo[k] = roundValueToNearestSIPrefix(latest - v[v.length - 8].value)
        carouselStats[k] = latest
    }

    const gridGap = '16px'

    return (
        <Container
            sx={{
                '& > *:not(:first-child, :last-child)': {
                    mb: itemGap,
                },
            }}
        >
            <Typography
                variant='h1'
                sx={{
                    ...headerMargins,
                }}
            >
                Library Stats
            </Typography>

            <Box
                sx={(theme) => ({
                    display: 'grid',
                    gridTemplateColumns: 'repeat(12, 1fr)',
                    columnGap: gridGap,
                    '&, & > .stat-graph': {
                        height: '653px',
                        [theme.breakpoints.up('tablet')]: {
                            height: '659px',
                        },
                        [theme.breakpoints.up('desktop')]: {
                            height: '659px',
                        },
                    },
                    '& .stat-graph': {
                        gridColumn: 'span 12',
                        [theme.breakpoints.up('tablet')]: {
                            gridColumn: 'span 8',
                        },
                        [theme.breakpoints.up('desktop')]: {
                            gridColumn: 'span 9',
                        },
                    },
                    '& .stat-card': {
                        display: 'none',
                        [theme.breakpoints.up('tablet')]: {
                            display: 'flex',
                            gridColumn: 'span 4',
                        },
                        [theme.breakpoints.up('desktop')]: {
                            gridColumn: 'span 3',
                        },
                    },
                })}
            >
                <StatTimeseriesLineChart
                    data={stats.users}
                    dataLabelSingular='user'
                    dataLabelPlural='users'
                    headingLabel='Overall Membership'
                    headingValue={latestStats.users.rawValue}
                />
                <StatCard
                    primary={{
                        value: `${latestStats.bytes.value}${latestStats.bytes.SIUnitPrefix}B`,
                        label: 'Total Data Stored',
                    }}
                    secondary={{
                        value: `${since7DaysAgo.bytes.value}${since7DaysAgo.bytes.SIUnitPrefix}B`,
                        label: 'added this week',
                    }}
                    deltaSign='+'
                    textColor={darkText}
                    backgroundColor={theme.palette.yellow.grade50}
                    alwaysShowSecondary={true}
                />
            </Box>
            <Box
                sx={(theme) => ({
                    height: '198px',
                    [theme.breakpoints.up('tablet')]: {
                        display: 'none',
                    },
                    '& .stats-carousel': {
                        height: '100%',
                    },
                })}
            >
                <StatsCarousel stats={carouselStats} />
            </Box>
            <Box
                sx={(theme) => ({
                    display: 'grid',
                    // TODO?: switch to gridAutoColumns
                    // to avoid switching b/w gap definitions
                    gridTemplateColumns: 'repeat(13, 1fr)',
                    columnGap: itemGap,
                    [theme.breakpoints.up('tablet')]: {
                        gap: itemGap,
                    },
                })}
            >
                <Box
                    sx={(theme) => ({
                        display: 'grid',
                        gridTemplateColumns: 'repeat(2, 1fr)',
                        gap: itemGap,
                        [theme.breakpoints.up('tablet')]: {
                            gridColumn: 'span 13',
                        },
                        [theme.breakpoints.up('desktop')]: {
                            gridColumn: 'span 6',
                        },
                        '& .stat-card': {
                            display: 'none',
                            [theme.breakpoints.up('tablet')]: {
                                display: 'flex',
                                height: '284px',
                            },
                            [theme.breakpoints.up('desktop')]: {
                                height: '290px',
                            },
                        },
                    })}
                >
                    <StatCard
                        primary={{
                            value: latestStats.files.rawValue.toLocaleString(),
                            label: 'Files',
                        }}
                        secondary={{
                            value: since7DaysAgo.files.rawValue.toLocaleString(),
                            label: 'added this week',
                        }}
                        deltaSign='+'
                        deltaColor='light'
                        textColor={lightText}
                        backgroundColor={theme.palette.teal.grade50}
                    />
                    <StatCard
                        primary={{
                            value: latestStats.samples.rawValue.toLocaleString(),
                            label: 'Samples',
                        }}
                        secondary={{
                            value: since7DaysAgo.samples.rawValue.toLocaleString(),
                            label: 'added this week',
                        }}
                        deltaSign='+'
                        textColor={darkText}
                        backgroundColor={theme.palette.green.grade75}
                    />
                    <StatCard
                        primary={{
                            value: latestStats.subjects.rawValue.toLocaleString(),
                            label: 'Subjects',
                        }}
                        secondary={{
                            value: since7DaysAgo.subjects.rawValue.toLocaleString(),
                            label: 'added this week',
                        }}
                        deltaSign='+'
                        deltaColor='light'
                        textColor={lightText}
                        backgroundColor={theme.palette.orange.grade50}
                    />
                    <StatCard
                        primary={{
                            value: latestStats.assays.rawValue.toLocaleString(),
                            label: 'Assays',
                        }}
                        secondary={{
                            value: since7DaysAgo.assays.rawValue.toLocaleString(),
                            label: 'added this week',
                        }}
                        deltaSign='+'
                        deltaColor='light'
                        textColor={lightText}
                        backgroundColor={theme.palette.blue.grade50}
                    />
                </Box>
                <Box
                    sx={(theme) => ({
                        gridColumn: 'span 13',
                        [theme.breakpoints.up('desktop')]: {
                            gridColumn: 'span 7',
                        },
                    })}
                >
                    <ThemeProjectBreakdownChart
                        data={themeProjectBreakdown}    
                    />
                </Box>
            </Box>
        </Container>
    )
}