'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import ButtonBase from '@mui/material/ButtonBase';
import { alpha } from '@mui/material/styles'
import { useTheme } from '@mui/material/styles';
import { Swiper, SwiperSlide, SwiperClass } from 'swiper/react';
import { A11y, Autoplay } from 'swiper/modules'

import StatCard from './stat-card'


export type Stats = Record<'bytes' | 'assays' | 'subjects' | 'files' | 'samples' | 'users', number>


export default function StatsCarousel({ stats }: { stats: Stats }) {
    const theme = useTheme()
    const darkText = theme.palette.ground.grade10
    const lightText = theme.palette.utilityHighlight.main
    const bgTransition = theme.transitions.create(
        ['background-color'],
        {
            duration: theme.transitions.duration.quint,
            easing: theme.transitions.easing.quint,
        },
    )
    const transformTransition = theme.transitions.create(
        ['transform'],
        {
            duration: theme.transitions.duration.quint,
            easing: theme.transitions.easing.quint,
        },
    )

    const [itemIndex, setItemIndex] = React.useState<number>(0)
    const [swiperRef, setSwiperRef] = React.useState<SwiperClass>()

    const items = [
        {
            label: 'Total Data Stored',
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
            label: 'Samples',
            value: stats.samples.toLocaleString(),
            textColor: darkText,
            backgroundColor: theme.palette.green.grade75,
        },

        {
            label: 'Users',
            value: stats.users.toLocaleString(),
            textColor: darkText,
            backgroundColor: theme.palette.magenta.grade75,
        },
    ]

    const handleClickCarouselIndexIndicator = (index: number) => {
        setItemIndex(index)
        swiperRef?.slideTo(index, undefined, false)
    }

    const simpleStatPaddingPx = 16

    return (
        <Box
            className='stats-carousel'
            sx={{
                display: 'flex',
                flexDirection: 'column',
                justifyContent: 'space-between',
                borderRadius: '30px',
                backgroundColor: items[itemIndex].backgroundColor,
                transition: bgTransition,
                '& > *:first-child': {
                    flexGrow: 1,
                },
            }}
        >
            <Box
                sx={{
                    '& .swiper, .swiper-wrapper, .swiper-slide': {
                        height: '100%',
                    },
                    '& .stat-card': {
                        // TODO: why can't this be just 100%?
                        height: `calc(100% - ${simpleStatPaddingPx * 2}px)`,
                        padding: `${simpleStatPaddingPx * 2}px)`
                    },
                    '& .swiper-wrapper': {
                        transition: transformTransition,
                    },
                }}
            >
                <Swiper
                    modules={[A11y, Autoplay]}
                    autoplay={{
                        delay: theme.transitions.duration.long,
                    }}
                    onInit={(swiper) => setSwiperRef(swiper)}
                    onSlideChange={(swiper) => setItemIndex(swiper.activeIndex)}
                >
                    {items.map((item) => {
                        return (
                            <SwiperSlide key={item.label}>
                                <StatCard
                                    primary={{
                                        value: item.value,
                                        label: item.label,
                                    }}
                                    textColor={item.textColor}
                                />
                            </SwiperSlide>
                        )
                    })}
                </Swiper>
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
                        key={index}
                        className={`carousel-index-indicator${isActive ? ' active' : ''}`}
                        onClick={() => handleClickCarouselIndexIndicator(index)}
                    />
                })}
            </Box>
        </Box >
    )
}