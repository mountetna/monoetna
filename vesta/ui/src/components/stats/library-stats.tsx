'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Container from '@mui/system/Container'
import Typography from '@mui/material/Typography';

import StatCard from './stat-card'
import StatsCarousel, { Stats as StatsProp } from './stats-carousel'
import StatGraph from './stat-graph'
import { Instance } from './types';
import theme, { headerMargins } from '@/theme';
import { SIValue, roundValueToNearestSIUnit } from '@/lib/utils/units';


export type StatsTimeseries = Record<keyof StatsProp, Instance<number>[]>

export default function LibraryStats({ stats }: { stats: StatsTimeseries }) {
    const gridGap = '16px'

    // @ts-ignore
    const latestStats: Record<keyof StatsProp, SIValue> = {}
    // @ts-ignore
    const since7DaysAgo: Record<keyof StatsProp, SIValue> = {}

    for (const [k, v] of Object.entries(stats)) {

        const latest = v[v.length - 1].value
        // @ts-ignore
        latestStats[k] = roundValueToNearestSIUnit(latest)
        // @ts-ignore
        since7DaysAgo[k] = roundValueToNearestSIUnit(latest - v[v.length - 8].value)
    }

    return (
        <Container>
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
                    '& .stat-graph': {
                        gridColumn: 'span 12',
                        [theme.breakpoints.up('tablet')]: {
                            gridColumn: 'span 8',
                        },
                        [theme.breakpoints.up('desktop')]: {
                            gridColumn: 'span 9',
                        },
                    },
                    '& .simple-stat': {
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
                <StatGraph />
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
                    textColor='utilityLowlight.main'
                    backgroundColor={theme.palette.yellow.grade50}
                />
            </Box>
            <Box
                sx={(theme) => ({
                    [theme.breakpoints.up('tablet')]: {
                        display: 'none',
                    },
                })}
            >
                <StatsCarousel stats={latestStats} />
            </Box>
            <Box
                sx={(theme) => ({
                    display: 'none',
                    [theme.breakpoints.up('tablet')]: {
                        display: 'grid',
                        gridTemplateColumns: 'repeat(2, 1fr)',
                        columnGap: '',
                    },
                })}
            >

            </Box>
        </Container>
    )
}