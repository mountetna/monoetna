'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import ButtonBase from '@mui/material/ButtonBase';
import { alpha } from '@mui/material/styles'
import { useTheme } from '@mui/material/styles';
import Carousel from 'react-multi-carousel'

import SimpleStat from '../stats/simple-stat'


export interface Stats {
    bytes: number
    assays: number
    subjects: number
    files: number
    users: number
}


export default function StatsCarousel({ stats }: { stats: Stats }) {
    const theme = useTheme()
    const darkText = theme.palette.ground.grade10
    const lightText = theme.palette.utilityHighlight.main

    const [itemIndex, setItemIndex] = React.useState<number>(0)
    const carouselRef = React.createRef<Carousel>()

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

    const handleClickCarouselIndexIndicator = (index: number) => {
        setItemIndex(index)
        carouselRef.current?.goToSlide(index, true)
    }

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

    return (
        <Box
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
                    '& .react-multi-carousel-list, .react-multi-carousel-track, .simple-stat': {
                        height: '100%',
                    },
                }}
            >
                <Carousel
                    ref={carouselRef}
                    responsive={{
                        global: {
                            breakpoint: { max: 3000, min: 0 },
                            items: 1,
                        }
                    }}
                    beforeChange={(nextSlide) => setItemIndex(nextSlide)}
                    swipeable={true}
                    draggable={true}
                    keyBoardControl={true}
                    ssr={true}
                    autoPlay={true}
                    autoPlaySpeed={2000}
                    customTransition={transformTransition}
                    arrows={false}
                    rewind={true}
                    rewindWithAnimation={true}
                >
                    {items.map((item) => {
                        return (
                            <React.Fragment key={item.label}>
                                <SimpleStat
                                    primary={{
                                        value: item.value,
                                        label: item.label,
                                    }}
                                    textColor={item.textColor}
                                />
                            </React.Fragment>
                        )
                    })}
                </Carousel>
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