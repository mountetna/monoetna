import * as React from 'react';
import { Metadata } from 'next';
import Box from '@mui/system/Box'

import Hero from '@/components/hero/hero';
import { getRandomItem } from '@/lib/utils/random';
import AboutCarousel from '@/components/about/about-carousel';
import LibraryStats from '@/components/stats/library-stats';
import { sectionMargins } from '@/theme';

import oscc1Fallback from '/public/images/hero/oscc1-fallback.png'
import xeniumFallback from '/public/images/hero/xenium-fallback.png'
import tonsilFallback from '/public/images/hero/tonsil-fallback.png'
import aboutImg1 from '/public/images/about-carousel/carousel-img-1.png'
import aboutImg2 from '/public/images/about-carousel/carousel-img-2.png'
import aboutImg3 from '/public/images/about-carousel/carousel-img-3.png'
import aboutImg4 from '/public/images/about-carousel/carousel-img-4.png'
import { Instance } from '@/components/stats/types';


export const metadata: Metadata = {
  title: 'UCSF Data Library',
}

const VIDEOS = [
  { videoSrc: '/videos/hero/oscc1.mp4', imageSrc: oscc1Fallback },
  { videoSrc: '/videos/hero/xenium.mp4', imageSrc: xeniumFallback },
  { videoSrc: '/videos/hero/tonsil.mp4', imageSrc: tonsilFallback },
]

const ABOUT_ITEMS = [
  {
    title: 'About the Library',
    header: 'About the Data Library',
    body: 'The UCSF Data Library project aims to capture, curate, and share biological data generated on campus — enabling data search, exploration, visualization and cross-project analyses through a series of applications in available in your web browser. The library engineering team in the Data Science CoLab is building these tools with close collaboration and input from the CoLabs, ImmunoX, and participating CoProject labs.',
    image: {
      src: aboutImg1,
      alt: 'Random Data Library iconography',
    },
  },
  {
    title: 'Contributing',
    header: 'Contributing to the library',
    body: 'Anyone at UCSF can contribute data from their project, and external collaborators can reach out to add their data from the library while sharing data alike. The library is here to be a source for open collaboration to enable new research.',
    link: {
      header: 'Working on something great?',
      blurb: 'We’re always looking for new data and projects to make more accessible.',
      href: '#',
      label: 'Get in touch',
    },
    image: {
      src: aboutImg2,
      alt: 'Biologist holding a test tube',
    },
  },
  {
    title: 'Available tools',
    header: 'What tools are available',
    body: 'The data library has many tools at your disposal, Janus for administration and permissions, midas for enim sem ut tincidunt vehicula vulputate donec etiam morbi. Ac nec facilisis sagittis aliquet. Felis laoreet sed rhoncus a quis odio dignissim est nisl. At est vulputate id etiam felis. Proin ac dapibus nec a pretium vel. Tincidunt feugiat dolor risus maecenas est est varius senectus scelerisque.',
    link: {
      href: '#',
      label: 'Get Access',
    },
    image: {
      src: aboutImg3,
      alt: 'Screenshot of Data Library tools',
    },
  },
  {
    title: 'Getting access',
    header: 'Getting access to the library',
    body: 'The Data Library requires a library card to access. Each project is unique, and must be added to your card separately. All projects require signing a data sharing agreement to access.',
    link: {
      href: '#',
      label: 'Get Access',
    },
    image: {
      src: aboutImg4,
      alt: 'A few Data Library library cards',
    },
  },
]

function getRandomArbitrary(min: number, max: number): number {
  return Math.floor(Math.random() * (max - min) + min)
}

const startDate = new Date(Date.now())
startDate.setDate(startDate.getDate() - 365)

const STATS = {
  bytes: Array(365).fill(null).map((_, i) => ({
    timestamp: (new Date(startDate.getDate() + i)).getTime(),
    value: getRandomArbitrary(1e12, 999e12),
  })),
  assays: Array(365).fill(null).map((_, i) => ({
    timestamp: (new Date(startDate.getDate() + i)).getTime(),
    value: getRandomArbitrary(1000, 400000),
  })),
  subjects: Array(365).fill(null).map((_, i) => ({
    timestamp: (new Date(startDate.getDate() + i)).getTime(),
    value: getRandomArbitrary(10, 999),
  })),
  files: Array(365).fill(null).map((_, i) => ({
    timestamp: (new Date(startDate.getDate() + i)).getTime(),
    value: getRandomArbitrary(1000, 999999),
  })),
  samples: Array(365).fill(null).map((_, i) => ({
    timestamp: (new Date(startDate.getDate() + i)).getTime(),
    value: getRandomArbitrary(100000, 999999),
  })),
  users: Array(365).fill(null).map((_, i) => ({
    timestamp: (new Date(startDate.getDate() + i)).getTime(),
    value: getRandomArbitrary(100, 1000),
  })),
}

async function getData() {
  return {
    heroVideo: getRandomItem(VIDEOS),
    // TODO: replace with real data
    stats: STATS,
    // TODO: replace with real data
    aboutItems: ABOUT_ITEMS,
  }
}

export default async function Home() {
  const data = await getData()

  const carouselStats = {} as Record<keyof typeof STATS, number>
  for (const [k, v] of (Object.entries(STATS) as [keyof typeof STATS, Instance<number>[]][])) {
    carouselStats[k] = v[v.length - 1].value
  }

  return (
    <React.Fragment>
      <Hero
        video={data.heroVideo}
        stats={carouselStats}
        scrollTargetId='about'
      />

      <Box sx={sectionMargins}>
        <Box id='about'>
          <AboutCarousel
            items={data.aboutItems}
          />
        </Box>

        <Box id='stats'>
          <LibraryStats
            stats={data.stats}
          />
        </Box>
      </Box>

    </React.Fragment>
  );
}
