'use client'

import * as React from 'react';
import Box from '@mui/system/Box'
import Container from '@mui/system/Container';
import Image from 'next/image';
import Typography from '@mui/material/Typography';
import ButtonBase from '@mui/material/ButtonBase';

import oscc1Fallback from '/public/images/hero/oscc1-fallback.png'
import xeniumFallback from '/public/images/hero/xenium-fallback.png'
import arrowUpRightLightSrc from '/public/images/icons/arrow-up-right-light.svg'

import { getRandomItem } from '@/lib/utils/random';


// TODO: fix nextjs server vs client src prop disagreement
const VIDEOS = [
    { videoSrc: '/videos/hero/oscc1.mp4', imageSrc: oscc1Fallback },
    { videoSrc: '/videos/hero/xenium.mp4', imageSrc: xeniumFallback },
]

export default function Hero() {
    const video = getRandomItem(VIDEOS)

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
                        [theme.breakpoints.up('desktop')]: {
                            gridTemplateColumns: 'repeat(14, 1fr)',
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
                            [theme.breakpoints.up('tablet')]: {
                                gridColumn: 'span 10',
                            },
                            [theme.breakpoints.up('desktop')]: {
                                gridColumn: 'span 7',
                            },
                        })}
                    >
                        <video
                            // poster={video.imageSrc.src}
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
                            />
                        </video>
                    </Box>
                    <Box
                        sx={(theme) => ({
                            gridColumn: 'span 12',
                            color: 'utilityHighlight.main',
                            [theme.breakpoints.up('tablet')]: {
                                gridColumn: 'span 10',
                            },
                            [theme.breakpoints.up('desktop')]: {
                                gridColumn: 'span 7',
                            },
                        })}
                    >
                        <Typography variant='h3'>
                            A library capturing, curating, and sharing biological data
                            generated on UCSF campus—enabling radically collaborative research.
                        </Typography>
                        <Box>
                            <ButtonBase
                            // sx={{
                            //     px: '16px',
                            //     py: '8px',
                            //     color: 'utilityHighlight.main',
                            //     backgroundColor: 'blue.grade50',
                            //     borderRadius: '10px',
                            // }}
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
                </Box>
            </Container>
        </Box>
    )
}