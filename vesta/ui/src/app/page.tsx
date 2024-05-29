import * as React from 'react';
import { Metadata } from 'next';
import Typography from '@mui/material/Typography';
import Box from '@mui/system/Box';
import NextLink from 'next/link';
import Container from '@mui/system/Container';

import Hero from '@/components/hero/Hero';
import { getRandomItem } from '@/lib/utils/random';
import oscc1Fallback from '/public/images/hero/oscc1-fallback.png'
import xeniumFallback from '/public/images/hero/xenium-fallback.png'


export const metadata: Metadata = {
  title: 'UCSF Data Library',
}

const VIDEOS = [
  { videoSrc: '/videos/hero/oscc1.mp4', imageSrc: oscc1Fallback.src },
  { videoSrc: '/videos/hero/xenium.mp4', imageSrc: xeniumFallback.src },
]

async function getData() {
  return {
    heroVideo: getRandomItem(VIDEOS),
    // TODO: replace with real data
    stats: {
      bytes: 376e12,
      assays: 348500,
      subjects: 221,
      files: 641202,
      users: 519,
    }
  }
}

export default async function Home() {
  const data = await getData()

  return (
    <React.Fragment>
      <Hero
        video={data.heroVideo}
        stats={data.stats}
      />
    </React.Fragment>
  );
}
