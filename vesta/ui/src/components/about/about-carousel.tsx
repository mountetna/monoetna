'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import { useTheme } from '@mui/material/styles';
import useMediaQuery from '@mui/material/useMediaQuery';
import ButtonBase from '@mui/material/ButtonBase';
import { Swiper, SwiperSlide, SwiperClass } from 'swiper/react';
import { A11y, EffectFade } from 'swiper/modules'

import AboutItem, { Link, ImageProps } from './about-item'
import Typography from '@mui/material/Typography';


interface AboutItemProps {
    title: string
    header: string
    body: string
    link?: Link
    image: ImageProps
}


export default function AboutCarousel({ items }: { items: AboutItemProps[] }) {
    const theme = useTheme()
    const isMobile = useMediaQuery(theme.breakpoints.between(
        theme.breakpoints.values.mobile,
        theme.breakpoints.values.tablet,
    ))

    const indexIndicatorFadeMaskWidthPx = 36

    const transition = theme.transitions.create(
        ['all'],
        {
            duration: theme.transitions.duration.quint,
            easing: theme.transitions.easing.quint,
        },
    )

    const [swiperRefMobile, setSwiperRefMobile] = React.useState<SwiperClass>()
    const [swiperRefTabletDesktop, setSwiperRefTabletDesktop] = React.useState<SwiperClass>()
    const [itemIndex, setItemIndex] = React.useState(0)
    const indexIndicatorsContainer = React.createRef<HTMLElement>()
    const indexIndicatorRefs = items.map(_ => React.createRef<HTMLElement>())

    const handleClickCarouselIndexIndicator = (index: number) => {
        setItemIndex(index)

        swiperRefMobile?.slideTo(index, undefined, false)
        swiperRefTabletDesktop?.slideTo(index, undefined, false)

        if (isMobile) {
            const container = indexIndicatorsContainer.current
            const indexIndicator = indexIndicatorRefs[index].current

            if (indexIndicator === null || container === null) { return }

            const indexIndicatorPosition = indexIndicator.offsetLeft

            container.scrollTo({
                left: indexIndicatorPosition - indexIndicatorFadeMaskWidthPx,
                behavior: 'smooth',
            })
        }
    }

    return (
        <Box
            sx={{
                '& .swiper-mobile': {
                    [theme.breakpoints.up('tablet')]: {
                        display: 'none',
                    },
                },
                '& .swiper-tablet-desktop': {
                    display: 'none',
                    [theme.breakpoints.up('tablet')]: {
                        display: 'block',
                    },
                },
                '& .swiper-wrapper': {
                    transition: transition
                },
                '& .swiper-slide': {
                    width: 'fit-content',
                },
            }}
        >
            <Box
                sx={{
                    position: 'relative',
                    mb: '24px',
                    mx: '8px',
                    [theme.breakpoints.up('tablet')]: {
                        mb: '37px',
                    },
                    [theme.breakpoints.up('desktop')]: {
                        mb: '37px',
                    },
                    '& .fade-mask': {
                        position: 'absolute',
                        zIndex: 1,
                        top: 0,
                        width: `${indexIndicatorFadeMaskWidthPx}px`,
                        height: '100%',
                        [theme.breakpoints.up('tablet')]: {
                            display: 'none',
                        },
                    },
                }}
            >
                <Box
                    className='fade-mask'
                    sx={{
                        left: 0,
                        background: `linear-gradient(to left, transparent, ${theme.palette.utilityHighlight.main})`
                    }}
                />
                <Box
                    ref={indexIndicatorsContainer}
                    sx={{
                        display: 'flex',
                        overflow: 'scroll',
                        textWrap: 'nowrap',
                        px: `${indexIndicatorFadeMaskWidthPx}px`,
                        [theme.breakpoints.up('tablet')]: {
                            justifyContent: 'center',
                            overflow: 'unset',
                            width: 'auto',
                        },
                        [theme.breakpoints.up('desktop')]: {
                        },
                        '& .carousel-index-indicator-container': {
                            '&:not(:last-child)': {
                                mr: '17px',
                                [theme.breakpoints.up('tablet')]: {
                                    mr: '60px',
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
                            transition: transition,
                            '&.active, &:hover': {
                                color: 'blue.grade50',
                            },
                            '&.active': {
                                borderColor: 'blue.grade50',
                            },
                        },
                    }}
                >
                    {items.map((item, index) => {
                        const isActive = index === itemIndex
                        return (
                            <Box
                                ref={indexIndicatorRefs[index]}
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
                <Box
                    className='fade-mask'
                    sx={{
                        right: 0,
                        background: `linear-gradient(to right, transparent, ${theme.palette.utilityHighlight.main})`
                    }}
                />
            </Box>
            <Swiper
                className='swiper-mobile'
                modules={[A11y, EffectFade]}
                allowTouchMove={false}
                autoHeight={true}
                onInit={(swiper) => setSwiperRefMobile(swiper)}
                onSlideChange={(swiper) => setItemIndex(swiper.activeIndex)}
                effect={'fade'}
                fadeEffect={{ crossFade: true, }}
                speed={theme.transitions.duration.quint}
                spaceBetween={0}
            >
                {items.map((item) => {
                    return (
                        <SwiperSlide key={item.header}>
                            <AboutItem
                                header={item.header}
                                body={item.body}
                                link={item.link}
                                image={item.image}
                            />
                        </SwiperSlide>
                    )
                })}
            </Swiper>
            <Box
                sx={{
                    position: 'relative',
                    '& .fade-mask': {
                        display: 'none',
                        position: 'absolute',
                        zIndex: 10,
                        top: 0,
                        width: '47px',
                        height: '100%',
                        [theme.breakpoints.up(theme.breakpoints.values.desktopLg + 47)]: {
                            display: 'block',
                        },
                        [theme.breakpoints.up(theme.breakpoints.values.desktopLg + 300)]: {
                            width: '220px'
                        },
                    },
                }}
            >
                <Box
                    className='fade-mask'
                    sx={{
                        left: 0,
                        background: `linear-gradient(to left, transparent, ${theme.palette.utilityHighlight.main})`
                    }}
                />
                <Swiper
                    className='swiper-tablet-desktop'
                    modules={[A11y]}
                    allowTouchMove={false}
                    autoHeight={true}
                    onInit={(swiper) => setSwiperRefTabletDesktop(swiper)}
                    onSlideChange={(swiper) => setItemIndex(swiper.activeIndex)}
                    centeredSlides={true}
                    slidesPerView={'auto'}
                    speed={theme.transitions.duration.quint}
                    breakpoints={{
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
                                    image={item.image}
                                />
                            </SwiperSlide>
                        )
                    })}
                </Swiper>
                <Box
                    className='fade-mask'
                    sx={{
                        right: 0,
                        background: `linear-gradient(to right, transparent, ${theme.palette.utilityHighlight.main})`
                    }}
                />
            </Box>
        </Box>
    )
}