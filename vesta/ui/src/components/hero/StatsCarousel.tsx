'use client'

import * as React from 'react'
import SwipeableViews from 'react-swipeable-views'
import { autoPlay, bindKeyboard } from 'react-swipeable-views-utils'

import SimpleStat from '../stats/SimpleStat'
import theme from '@/theme'


const EnhancedSwipeableViews = autoPlay(bindKeyboard(SwipeableViews))

export interface Stats {
    bytes: number
    assays: number
    subjects: number
    files: number
    users: number
}

export default function StatsCarousel({ stats }: { stats: Stats }) {
    return (
        <EnhancedSwipeableViews
            enableMouseEvents
            autoPlay={true}
            interval={theme.transitions.duration.long}
        >
            <SimpleStat
                primary={{
                    value: `${Math.round(stats.bytes / 1e12).toLocaleString()}TB`,
                    label: 'Total Data',
                }}
            />
            <SimpleStat
                primary={{
                    value: stats.assays.toLocaleString(),
                    label: 'Assays',
                }}
            />
            <SimpleStat
                primary={{
                    value: stats.subjects.toLocaleString(),
                    label: 'Subjects',
                }}
            />
            <SimpleStat
                primary={{
                    value: stats.files.toLocaleString(),
                    label: 'Files',
                }}
            />
            <SimpleStat
                primary={{
                    value: stats.users.toLocaleString(),
                    label: 'Users',
                }}
            />
        </EnhancedSwipeableViews>
    )
}