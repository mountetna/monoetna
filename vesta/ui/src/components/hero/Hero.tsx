'use client'

import * as React from 'react';
import Box from '@mui/system/Box'
import Container from '@mui/system/Container';
import Image from 'next/image';
import Typography from '@mui/material/Typography';
import ButtonBase from '@mui/material/ButtonBase';
import MUILink from '@mui/material/Link';
import Link from 'next/link'

import arrowUpRightLightSrc from '/public/images/icons/arrow-up-right-light.svg'
import StatsCarousel, { Stats } from './StatsCarousel';


export interface Video {
    videoSrc: string
    imageSrc: string
}

export default function Hero({ video, stats }: { video: Video, stats: Stats }) {
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
                            poster={video.imageSrc}
                            width='100%'
                            height='auto'
                            playsInline
                            autoPlay
                            loop
                            muted
                        >
                            <source
                                src={video.videoSrc}
                                type='video/mp4'
                            />
                            <Image
                                src={video.imageSrc}
                                alt='Picture of Spatial Transcriptomics'
                                width={1080}
                                height={1080}
                            />
                        </video>
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
                            A library capturing, curating, and sharing biological data
                            generated on UCSF campus—enabling radically collaborative research.
                        </Typography>
                        <Box>
                            <ButtonBase>
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
                                        src={arrowUpRightLightSrc}
                                        alt='sda'
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
                                    height: '211px',
                                },
                                [theme.breakpoints.up('desktop')]: {
                                    height: '211px',
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
                            href="#"
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
                                Get Access
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
                                        src={arrowUpRightLightSrc}
                                        alt='sda'
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