import * as React from 'react';
import Box from '@mui/system/Box'
import Container from '@mui/system/Container';
import Image from 'next/image';

// import fallbackImg from 'public/images/hero/oscc1_fallback.png'


export default function Hero() {
    return (
        <Box
            sx={{
                backgroundColor: 'black',
            }}
        >
            <Container maxWidth='desktopLg'>
                <video
                    width={500}
                    height={500}
                    preload='none'
                    playsInline
                    autoPlay
                >
                    <source
                        src='/videos/hero/oscc1.mp4'
                        type='video/mp4'
                    />
                    <Image
                        src='/images/hero/oscc1_fallback.png'
                        alt='Picture of Spatial Transcriptomics'
                        width={500}
                        height={500}
                    />
                </video>
            </Container>
        </Box>
    )
}