'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import { useTheme } from '@mui/material/styles';
import ButtonBase from '@mui/material/ButtonBase';
import { Swiper, SwiperSlide, SwiperClass } from 'swiper/react';
import { A11y } from 'swiper/modules'

import AboutItem, { Link } from './about-item'
import Typography from '@mui/material/Typography';


interface AboutItemProps {
    title: string
    header: string
    body: string
    link?: Link
    imageSrc: string
}

export default function AboutCarousel({ items }: { items: AboutItemProps[] }) {
    const theme = useTheme()
    const navTransition = theme.transitions.create(
        ['all'],
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

    const [swiperRef, setSwiperRef] = React.useState<SwiperClass>()
    const [itemIndex, setItemIndex] = React.useState(0)

    const handleClickCarouselIndexIndicator = (index: number) => {
        setItemIndex(index)
        swiperRef?.slideTo(index)
    }

    return (
        <Box
            sx={{
                '& .swiper-slide': {
                    width: 'fit-content',
                },
            }}
        >
            <Box
                sx={{
                    display: 'flex',
                    justifyContent: 'center',
                    // width: 'fit-content',
                    overflow: 'scroll',
                    textWrap: 'nowrap',
                    mb: '24px',
                    [theme.breakpoints.up('tablet')]: {
                        overflow: 'unset',
                        width: 'auto',
                        mb: '37px',
                    },
                    [theme.breakpoints.up('desktop')]: {
                        mr: '37px',
                    },
                    '& .carousel-index-indicator-container': {
                        '&:not(:last-child)': {
                            mr: '17px',
                            [theme.breakpoints.up('tablet')]: {
                                mr: '71px',
                            },
                            [theme.breakpoints.up('desktop')]: {
                                mr: '71px',
                            },
                        },
                    },
                    '& .carousel-index-indicator': {
                        color: 'ground.grade10',
                        pb: '8px',
                        borderBottom: '3px solid transparent',
                        transition: theme.transitions.create(
                            ['all'],
                            {
                                duration: theme.transitions.duration.ease,
                                easing: theme.transitions.easing.ease,
                            }
                        ),
                        '&.active, &:hover': {
                            color: 'blue.grade50',
                        },
                        '&.active': {
                            borderColor: 'blue.grade50',
                        },
                    }
                }}
            >
                {items.map((item, index) => {
                    const isActive = index === itemIndex
                    return (
                        <Box
                            key={index}
                            className='carousel-index-indicator-container'
                        >
                            <ButtonBase
                                tabIndex={0}
                                className={`carousel-index-indicator${isActive ? ' active' : ''}`}
                                onClick={() => handleClickCarouselIndexIndicator(index)}
                            >
                                <Typography variant='h6'>
                                    {item.title}
                                </Typography>
                            </ButtonBase>
                        </Box>
                    )
                })}
            </Box>
            <Swiper
                modules={[A11y]}
                allowTouchMove={false}
                onInit={(swiper) => setSwiperRef(swiper)}
                centeredSlides={true}
                slidesPerView={'auto'}
                breakpoints={{
                    [theme.breakpoints.values.mobile]: {
                        spaceBetween: 0,
                    },
                    [theme.breakpoints.values.tablet]: {
                        spaceBetween: 0,
                    },
                    [theme.breakpoints.values.desktop]: {
                        spaceBetween: 47,
                    },
                }}
            >
                {items.map((item) => {
                    return (
                        <SwiperSlide key={item.header}>
                            <AboutItem
                                header={item.header}
                                body={item.body}
                                link={item.link}
                                imageSrc={item.imageSrc}
                            />
                        </SwiperSlide>
                    )
                })}
            </Swiper>
        </Box>
    )
}