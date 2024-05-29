'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import ButtonBase from '@mui/material/ButtonBase';
import { alpha } from '@mui/material/styles'
import SwipeableViews from 'react-swipeable-views'
import { autoPlay, bindKeyboard, virtualize, SlideRenderProps } from 'react-swipeable-views-utils'
import { mod } from 'react-swipeable-views-core';

import SimpleStat from '../stats/SimpleStat'
import theme from '@/theme'


const EnhancedSwipeableViews = bindKeyboard(autoPlay(SwipeableViews))

export interface Stats {
    bytes: number
    assays: number
    subjects: number
    files: number
    users: number
}

const darkText = theme.palette.ground.grade10
const lightText = theme.palette.utilityHighlight.main

export default function StatsCarousel({ stats }: { stats: Stats }) {
    const [itemIndex, setItemIndex] = React.useState<number>(0)
    const [carouselIndex, setCarouselIndex] = React.useState<number>(0)
    const carouselContainerRef = React.createRef<HTMLElement>()
    const [carouselContainerHeight, setCarouselContainerHeight] = React.useState<number>()

    const items = [
        {
            label: 'Total Data',
            value: `${Math.round(stats.bytes / 1e12).toLocaleString()}TB`,
            textColor: darkText,
            backgroundColor: theme.palette.yellow.grade50,
        },
        {
            label: 'Assays',
            value: stats.assays.toLocaleString(),
            textColor: lightText,
            backgroundColor: theme.palette.blue.grade50,
        },
        {
            label: 'Subjects',
            value: stats.subjects.toLocaleString(),
            textColor: lightText,
            backgroundColor: theme.palette.orange.grade50,
        },
        {
            label: 'Files',
            value: stats.files.toLocaleString(),
            textColor: lightText,
            backgroundColor: theme.palette.teal.grade50,
        },
        {
            label: 'Users',
            value: stats.users.toLocaleString(),
            textColor: darkText,
            backgroundColor: theme.palette.magenta.grade75,
        },
    ]

    const containerRef = React.createRef<HTMLElement>()

    const handleChangeIndex = (indexNext: number, indexCurrent: number) => {
        if (containerRef.current === null) { return }

        setCarouselIndex(indexNext)
        setItemIndex(mod(indexNext, items.length))
    }

    // use with virtualize
    const slideRenderer = (params: SlideRenderProps) => {
        const { index, key } = params
        const itemIndex = mod(index, items.length)
        const item = items[itemIndex]

        return (
            <Box key={key}>
                <SimpleStat
                    primary={{
                        value: item.value,
                        label: item.label,
                    }}
                    textColor={item.textColor}
                />
            </Box>
        )
    }

    const bgTransition = theme.transitions.create(
        ['background-color'],
        {
            duration: theme.transitions.duration.quint,
            easing: theme.transitions.easing.quint,
        },
    )

    const handleClickCarouselIndexIndicator = (index: number) => {
        setCarouselIndex(index)
        setItemIndex(index)
    }

    React.useEffect(() => {
        setCarouselContainerHeight(carouselContainerRef.current?.offsetHeight)
    }, [carouselContainerRef.current])

    return (
        <Box
            ref={containerRef}
            sx={(theme) => ({
                display: 'flex',
                flexDirection: 'column',
                justifyContent: 'space-between',
                borderRadius: '30px',
                backgroundColor: items[itemIndex].backgroundColor,
                transition: bgTransition,
                '& > *:first-child': {
                    flexGrow: 1,
                },
            })}
        >
            <Box ref={carouselContainerRef}>
                <EnhancedSwipeableViews
                    enableMouseEvents
                    springConfig={{
                        duration: `${theme.transitions.duration.quint / 1000}s`,
                        easeFunction: theme.transitions.easing.quint,
                        delay: '0s',
                    }}
                    interval={theme.transitions.duration.long}
                    index={carouselIndex}
                    onChangeIndex={handleChangeIndex}
                >
                    {items.map((item) => {
                        return <SimpleStat
                            primary={{
                                value: item.value,
                                label: item.label,
                            }}
                            textColor={item.textColor}
                            heightPx={carouselContainerHeight}
                        />
                    })}
                </EnhancedSwipeableViews>
            </Box>
            <Box
                sx={{
                    display: 'inline-flex',
                    width: '100%',
                    justifyContent: 'center',
                    pb: '8px',
                    '& .carousel-index-indicator': {
                        display: 'inline-block',
                        width: '15px',
                        height: '15px',
                        borderRadius: '50%',
                        transition: bgTransition,
                        backgroundColor: alpha(theme.palette.ground.grade10, 0.4),
                        '&.active': {
                            backgroundColor: theme.palette.ground.grade10,
                        },
                        '&:not(:last-child)': {
                            mr: '11px',
                        },
                    },
                }}
            >
                {items.map((_, index) => {
                    const isActive = index === itemIndex
                    return <ButtonBase
                        className={`carousel-index-indicator${isActive ? ' active' : ''}`}
                        onClick={() => handleClickCarouselIndexIndicator(index)}
                    />
                })}
            </Box>
        </Box>
    )
}