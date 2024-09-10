'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import ButtonBase from '@mui/material/ButtonBase';
import { useTheme, alpha } from '@mui/material/styles';
import { Swiper, SwiperSlide, SwiperClass } from 'swiper/react';
import { A11y, Autoplay } from 'swiper/modules'
import Image from 'next/image';
import { Fade } from '@mui/material';

import StatCard from './stat-card'

import arrowCounterClockwiseLightIcon from '/public/images/icons/arrow-counter-clockwise-light.svg'


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
    const opacityTransition = theme.transitions.create(
        'opacity',
        {
            easing: theme.transitions.easing.ease,
            duration: theme.transitions.duration.ease,
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

    const handleClickRestart = () => {
        handleClickCarouselIndexIndicator(0)
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
                position: 'relative',
                overflow: 'hidden',
                '& > *:nth-child(1)': {
                    flexGrow: 1,
                },
            }}
        >
            <Box
                sx={{
                    position: 'relative',
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
                        stopOnLastSlide: true,
                    }}
                    onInit={(swiper) => setSwiperRef(swiper)}
                    onSlideChange={(swiper) => setItemIndex(swiper.activeIndex)}
                >
                    {items.map((item, idx) => {
                        return (
                            <SwiperSlide key={item.label}>
                                <Box
                                    sx={{
                                        position: 'relative',
                                        height: '100%',
                                    }}
                                >
                                    <StatCard
                                        primary={{
                                            value: item.value,
                                            label: item.label,
                                        }}
                                        textColor={item.textColor}
                                    />
                                </Box>
                            </SwiperSlide>
                        )
                    })}
                </Swiper>
            </Box>

            <Box
                sx={{
                    position: 'relative',
                    zIndex: 3,
                    display: 'inline-flex',
                    width: '100%',
                    justifyContent: 'center',
                    pb: '8px',
                    gap: '11px',
                    '& .carousel-index-indicator': {
                        display: 'inline-block',
                        width: '15px',
                        height: '15px',
                        borderRadius: '50%',
                        transition: bgTransition,
                        bgcolor: alpha(theme.palette.ground.grade10, 0.4),
                        '&.active': {
                            bgcolor: theme.palette.ground.grade10,
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

            {itemIndex === items.length - 1 && <Fade
                in={true}
                appear={true}
                easing={theme.transitions.easing.quint}
                timeout={theme.transitions.duration.quint}
            >
                <Box>
                    <Box
                        sx={{
                            position: 'absolute',
                            top: 0,
                            left: 0,
                            width: '100%',
                            height: '100%',
                            zIndex: 2,
                            '&:hover': {
                                '& > *:nth-child(1)': {
                                    opacity: 0.5,
                                },
                                '& > *:nth-child(2)': {
                                    opacity: 1,
                                },
                            },
                        }}
                    >
                        <Box
                            sx={{
                                position: 'absolute',
                                top: 0,
                                left: 0,
                                width: '100%',
                                height: '100%',
                                bgcolor: theme.palette.utilityLowlight.main,
                                opacity: 0.3,
                                transition: opacityTransition,
                            }}
                        />

                        <ButtonBase
                            onClick={handleClickRestart}
                            disabled={itemIndex !== items.length - 1}
                            aria-label='Restart'
                            sx={{
                                width: '100%',
                                height: '100%',
                                display: 'flex',
                                alignItems: 'center',
                                justifyContent: 'center',
                                opacity: 0.7,
                                transition: opacityTransition,
                            }}
                        >
                            <Image
                                src={arrowCounterClockwiseLightIcon}
                                alt='Arrow counter clockwise'
                                height={50}
                            />
                        </ButtonBase>
                    </Box>
                </Box>
            </Fade>}
        </Box>
    )
}