'use client'

import * as React from 'react';
import Box from '@mui/system/Box'
import Container from '@mui/system/Container';
import Image, { StaticImageData } from 'next/image';
import Typography from '@mui/material/Typography';
import ButtonBase from '@mui/material/ButtonBase';
import MUILink from '@mui/material/Link';
import { Fade, useTheme } from '@mui/material';
import Link from 'next/link'

import StatsCarousel, { Stats } from '@/components/stats/stats-carousel';
import { scrollTo } from '@/lib/utils/scroll';
import { useBreakpoint } from '@/lib/utils/responsive';
import { useUser } from '../user/context';

import arrowUpRightLightIcon from '/public/images/icons/arrow-up-right-light.svg'
import arrowRightLightIcon from '/public/images/icons/arrow-right-light.svg'
import arrowCounterClockwiseLightIcon from '/public/images/icons/arrow-counter-clockwise-light.svg'


export interface Video {
    videoSrc: string
    imageSrc: StaticImageData
}

export default function Hero({
    videos,
    initVideoIdx,
    stats,
    scrollTargetId,
    accessUrl,
}: {
    videos: Video[],
    initVideoIdx: number,
    stats: Stats,
    scrollTargetId: string,
    accessUrl: string,
}) {
    const theme = useTheme()
    const breakpoint = useBreakpoint()

    const user = useUser()

    const handleClickScrollToExplore = () => {
        const el = document.getElementById(scrollTargetId)

        if (el) {
            scrollTo({ top: el.offsetTop }, breakpoint)
        }
    }

    const [videoIdx, setVideoIdx] = React.useState(initVideoIdx)
    const video = videos[videoIdx]
    const videoRef = React.useRef<HTMLVideoElement>(null)

    const [videoPaused, setVideoPaused] = React.useState(false)
    const [videoEnded, setVideoEnded] = React.useState(false)

    const handleClickPlayPauseReset = async () => {
        const videoEl = videoRef.current
        if (!videoEl) return

        if (videoEnded) {
            videoEl.currentTime = 0
            await videoEl.play()
        } else if (videoPaused) {
            await videoEl.play()
        } else {
            videoEl.pause()
        }
    }

    const handleClickNext = async () => {
        const videoEl = videoRef.current
        if (!videoEl) return

        setVideoIdx(idx => (idx + 1) % videos.length)
        videoEl.load()
        await videoEl.play()
    }

    const videoControlButtonSizePx = 50

    const opacityTransition = theme.transitions.create(
        'opacity',
        {
            easing: theme.transitions.easing.ease,
            duration: theme.transitions.duration.ease,
        },
    )

    return (
        <Box
            sx={{
                backgroundColor: 'utilityLowlight.main',
            }}
        >
            <Container>
                <Box
                    sx={(theme) => ({
                        display: 'grid',
                        gridTemplateColumns: 'repeat(12, 1fr)',
                        gridTemplateRows: 'auto',
                        pt: '26px',
                        pb: '8px',
                        gap: '26px 0',
                        [theme.breakpoints.up('tablet')]: {
                            pb: '16px',
                        },
                        [theme.breakpoints.up('desktop')]: {
                            gridTemplateColumns: 'repeat(14, 1fr)',
                            alignItems: 'center',
                        },
                    })}
                >
                    <Box
                        sx={(theme) => ({
                            position: 'relative',
                            '& img': {
                                width: '100%',
                                height: 'auto',
                            },
                            gridColumn: 'span 12',
                            px: '31px',
                            [theme.breakpoints.up('tablet')]: {
                                gridColumn: 'span 10',
                                px: '0',
                            },
                            [theme.breakpoints.up('desktop')]: {
                                gridColumn: 'span 7',
                            },
                        })}
                    >
                        <video
                            ref={videoRef}
                            // poster={video.imageSrc.src}
                            playsInline
                            autoPlay
                            muted
                            onPlay={() => {
                                setVideoPaused(false)
                                setVideoEnded(false)
                            }}
                            onPause={() => setVideoPaused(true)}
                            onEnded={() => setVideoEnded(true)}
                            style={{
                                width: '100%',
                                height: 'auto',
                                aspectRatio: '1/1',
                                objectFit: 'cover',
                                borderRadius: '50%',
                            }}
                        >
                            <source
                                src={video.videoSrc}
                                type='video/mp4'
                            />

                            <Image
                                src={video.imageSrc}
                                alt='Picture of Spatial Transcriptomics'
                            />
                        </video>

                        {/* controls */}
                        <Box
                            onClick={() => videoEnded && handleClickPlayPauseReset()}
                            sx={{
                                position: 'absolute',
                                left: 0,
                                top: 0,
                                width: '100%',
                                height: '100%',
                                display: 'flex',
                                flexDirection: 'row',
                                alignItems: 'center',
                                justifyContent: 'center',
                                gap: '50px',
                                transition: opacityTransition,
                                opacity: 0,
                                '&:hover': {
                                    opacity: videoEnded ? 1 : 0,
                                    cursor: videoEnded ? 'pointer' : 'default',
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
                                    opacity: 0.7,
                                    bgcolor: theme.palette.utilityLowlight.main,
                                }}
                            />

                            {/* restart */}
                            <ButtonBase
                                onClick={handleClickPlayPauseReset}
                                disabled={!videoEnded}
                                aria-label='Restart video'
                                sx={{
                                    height: `${videoControlButtonSizePx}px`,
                                    aspectRatio: '1/1',
                                }}
                            >
                                <Image
                                    src={arrowCounterClockwiseLightIcon}
                                    alt='Arrow counter clockwise'
                                />
                            </ButtonBase>

                            {/* next */}
                            <ButtonBase
                                onClick={handleClickNext}
                                disabled={!videoEnded}
                                aria-label='Next video'
                                sx={{
                                    height: `${videoControlButtonSizePx}px`,
                                    aspectRatio: '1/1',
                                }}
                            >
                                <Image
                                    src={arrowRightLightIcon}
                                    alt='Arrow right'
                                />
                            </ButtonBase>
                        </Box>
                    </Box>
                    <Box
                        sx={(theme) => ({
                            gridColumn: 'span 12',
                            color: 'utilityHighlight.main',
                            py: '8px',
                            [theme.breakpoints.up('tablet')]: {
                                gridColumn: 'span 10',
                                px: '24px',
                                py: '24px',
                            },
                            [theme.breakpoints.up('desktop')]: {
                                gridColumn: 'span 7',
                            },
                        })}
                    >
                        <Typography
                            variant='h3'
                            sx={(theme) => ({
                                pb: '18px',
                            })}
                        >
                            A library capturing, curating, and sharing biological data — creating a foundation for transformative collaborative research.
                        </Typography>
                        <Box>
                            <ButtonBase
                                onClick={handleClickScrollToExplore}
                            >
                                <Box
                                    component='span'
                                    sx={{
                                        display: 'flex',
                                        justifyContent: 'center',
                                        alignItems: 'center',
                                        borderRadius: '50%',
                                        width: '52px',
                                        height: '52px',
                                        backgroundColor: 'ground.grade25',
                                        mr: '12px',
                                        '& img': {
                                            width: '32px',
                                            height: 'auto',
                                            transform: 'rotate(135deg)',
                                        }
                                    }}
                                >
                                    <Image
                                        src={arrowUpRightLightIcon}
                                        alt='Arrow pointing down'
                                    />
                                </Box>
                                <Typography variant='pMedium' color='ground.grade50'>
                                    Scroll to explore
                                </Typography>
                            </ButtonBase>
                        </Box>
                    </Box>
                    <Box
                        sx={(theme) => ({
                            gridColumn: 'span 12',
                            display: 'grid',
                            gridTemplateColumns: 'repeat(12, 1fr)',
                            gap: '8px',
                            alignContent: 'stretch',
                            color: 'utilityHighlight.main',
                            [theme.breakpoints.up('desktop')]: {
                                gridColumn: '8 / span 7',
                            },
                            '& > *': {
                                height: '211px',
                                [theme.breakpoints.up('tablet')]: {
                                    height: '203px',
                                },
                                [theme.breakpoints.up('desktop')]: {
                                    height: '198px',
                                },
                            },
                        })}
                    >
                        <Box
                            sx={(theme) => ({
                                gridColumn: 'span 12',
                                [theme.breakpoints.up('tablet')]: {
                                    gridColumn: 'span 6',
                                },
                                '& > *': {
                                    height: '100%',
                                },
                            })}
                        >
                            <StatsCarousel stats={stats} />
                        </Box>
                        <MUILink
                            href={accessUrl}
                            tabIndex={0}
                            component={Link}
                            underline='none'
                            sx={(theme) => ({
                                gridColumn: 'span 12',
                                display: 'flex',
                                flexDirection: 'column',
                                justifyContent: 'space-between',
                                p: '16px',
                                borderRadius: '30px',
                                backgroundColor: 'ground.grade25',
                                [theme.breakpoints.up('tablet')]: {
                                    gridColumn: 'span 6',
                                },
                            })}
                        >
                            <Typography
                                variant='h5'
                                color='utilityHighlight.main'
                            >
                                {user ? 'Go to Dashboard' : 'Get Access'}
                            </Typography>
                            <Box
                                sx={{
                                    display: 'flex',
                                    justifyContent: 'flex-end',
                                }}
                            >
                                <Box
                                    component='span'
                                    sx={{
                                        display: 'flex',
                                        justifyContent: 'center',
                                        alignItems: 'center',
                                        borderRadius: '50%',
                                        width: '64px',
                                        height: '64px',
                                        color: 'ground.grade75',
                                        backgroundColor: 'ground.grade10',
                                        '& img': {
                                            width: '39px',
                                            height: 'auto',
                                        }
                                    }}
                                >
                                    <Image
                                        src={arrowUpRightLightIcon}
                                        alt='Arrow pointing up-right'
                                    />
                                </Box>
                            </Box>
                        </MUILink>
                    </Box>
                </Box>
            </Container>
        </Box>
    )
}