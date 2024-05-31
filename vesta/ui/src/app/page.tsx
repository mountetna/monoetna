import * as React from 'react';
import { Metadata } from 'next';
import Box from '@mui/system/Box'

import Hero from '@/components/hero/hero';
import { getRandomItem } from '@/lib/utils/random';
import AboutCarousel from '@/components/about/about-carousel';
import { spacing } from '@/theme';

import oscc1Fallback from '/public/images/hero/oscc1-fallback.png'
import xeniumFallback from '/public/images/hero/xenium-fallback.png'
import tonsilFallback from '/public/images/hero/tonsil-fallback.png'
import aboutImg1 from '/public/images/about-carousel/carousel-img-1.png'
import aboutImg2 from '/public/images/about-carousel/carousel-img-2.png'
import aboutImg3 from '/public/images/about-carousel/carousel-img-3.png'
import aboutImg4 from '/public/images/about-carousel/carousel-img-4.png'


export const metadata: Metadata = {
  title: 'UCSF Data Library',
}

const VIDEOS = [
  { videoSrc: '/videos/hero/oscc1.mp4', imageSrc: oscc1Fallback.src },
  { videoSrc: '/videos/hero/xenium.mp4', imageSrc: xeniumFallback.src },
  { videoSrc: '/videos/hero/tonsil.mp4', imageSrc: tonsilFallback.src },
]

const ABOUT_ITEMS = [
  {
    title: 'About the Library',
    header: 'About the Data Library',
    body: 'The UCSF Data Library project aims to capture, curate, and share biological data generated on campus — enabling data search, exploration, visualization and cross-project analyses through a series of applications in available in your web browser. The library engineering team in the Data Science CoLab is building these tools with close collaboration and input from the CoLabs, ImmunoX, and participating CoProject labs.',
    imageSrc: aboutImg1.src,
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
    imageSrc: aboutImg2.src,
  },
  {
    title: 'Available tools',
    header: 'What tools are available',
    body: 'The data library has many tools at your disposal, Janus for administration and permissions, midas for enim sem ut tincidunt vehicula vulputate donec etiam morbi. Ac nec facilisis sagittis aliquet. Felis laoreet sed rhoncus a quis odio dignissim est nisl. At est vulputate id etiam felis. Proin ac dapibus nec a pretium vel. Tincidunt feugiat dolor risus maecenas est est varius senectus scelerisque.',
    link: {
      href: '#',
      label: 'Get Access',
    },
    imageSrc: aboutImg3.src,
  },
  {
    title: 'Getting access',
    header: 'Getting access to the library',
    body: 'The Data Library requires a library card to access. Each project is unique, and must be added to your card separately. All projects require signing a data sharing agreement to access.',
    link: {
      href: '#',
      label: 'Get Access',
    },
    imageSrc: aboutImg4.src,
  },
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
    },
    // TODO: replace with real data
    aboutItems: ABOUT_ITEMS,
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
      <Box sx={spacing}>
        <AboutCarousel items={data.aboutItems} />
      </Box>
    </React.Fragment>
  );
}
