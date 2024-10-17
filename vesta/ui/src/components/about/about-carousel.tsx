'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import { useTheme } from '@mui/material/styles';
import ButtonBase from '@mui/material/ButtonBase';
import { Swiper, SwiperSlide, SwiperClass } from 'swiper/react';
import { A11y, EffectFade } from 'swiper/modules'
import { usePathname, useRouter, useSearchParams } from 'next/navigation';
import Typography from '@mui/material/Typography';
import { useBreakpoint } from '@/lib/utils/responsive';

import AboutItem, { Link, ImageProps } from './about-item'
import { parseSearchParams, toSearchParamsString } from '@/lib/utils/uri';
import { ABOUT_SERACH_PARAMS_KEY, AboutSearchParamsState } from './models';


interface AboutItems {
    title: string
    header: string
    body: string
    link?: Link
    image: ImageProps
}

interface Props {
    items: AboutItems[]
}


export function _AboutCarousel({ items }: Props) {
    // Manage search params sync
    const router = useRouter()
    const pathname = usePathname()
    const searchParams = useSearchParams()

    React.useEffect(() => {
        const parsedSearchParams = parseSearchParams(searchParams)

        if (ABOUT_SERACH_PARAMS_KEY in parsedSearchParams) {
            const state: AboutSearchParamsState = parsedSearchParams[ABOUT_SERACH_PARAMS_KEY]

            if (state.index !== undefined) {
                handleClickCarouselIndexIndicator(state.index)
            }
        }
    }, [searchParams])

    const theme = useTheme()
    const isMobile = useBreakpoint() === 'mobile'

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
        updateUrl(index)

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

    const updateUrl = (itemIndex: number) => {
        const aboutState: AboutSearchParamsState = {
            index: itemIndex,
        }

        // push to router
        router.push(pathname + '?' + toSearchParamsString({ [ABOUT_SERACH_PARAMS_KEY]: aboutState }) + window.location.hash, { scroll: false })
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
                            '&.active': {
                                color: 'blue.grade50',
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
                initialSlide={itemIndex}
                className='swiper-mobile'
                modules={[A11y, EffectFade]}
                allowTouchMove={false}
                autoHeight={true}
                onInit={(swiper) => setSwiperRefMobile(swiper)}
                onSlideChange={(swiper) => setItemIndex(swiper.realIndex)}
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
                            width: '170px'
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
                    initialSlide={itemIndex}
                    className='swiper-tablet-desktop'
                    modules={[A11y]}
                    allowTouchMove={false}
                    autoHeight={true}
                    onInit={(swiper) => setSwiperRefTabletDesktop(swiper)}
                    onSlideChange={(swiper) => setItemIndex(swiper.realIndex)}
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


// Needed for https://nextjs.org/docs/messages/missing-suspense-with-csr-bailout
export default function AboutCarousel(props: Props) {
    return (
        <React.Suspense fallback={null}>
            <_AboutCarousel
                {...props}
            />
        </React.Suspense>
    )
}