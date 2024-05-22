'use client'

import * as React from 'react';
import Box from '@mui/system/Box'
import Container from '@mui/system/Container';
import Image from 'next/image';
import Typography from '@mui/material/Typography';
import ButtonBase from '@mui/material/ButtonBase';

import oscc1_fallback from '/public/images/hero/oscc1_fallback.png'

import { getRandomItem } from '@/lib/utils/random';


const VIDEOS = [
    '/videos/hero/oscc1.mp4',
    '/videos/hero/xenium.mp4',
]

export default function Hero() {


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
                            width='100%'
                            height='auto'
                            playsInline
                            autoPlay
                            loop
                            muted
                        >
                            <source
                                src={getRandomItem(VIDEOS)}
                                type='video/mp4'
                            />
                            <Image
                                src={oscc1_fallback}
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
                                //     typography: 'pBodyMediumWt',
                                // }}
                            >
                                <Box sx={{ backgroundColor: 'white' }} component='span'>AR</Box>
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