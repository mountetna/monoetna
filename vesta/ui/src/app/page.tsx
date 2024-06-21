import * as React from 'react';
import { Metadata } from 'next';
import Box from '@mui/system/Box'
import { faker } from '@faker-js/faker'
import _ from 'lodash'

import Hero from '@/components/hero/hero';
import { getRandomItem } from '@/lib/utils/random';
import AboutCarousel from '@/components/about/about-carousel';
import LibraryStats from '@/components/stats/library-stats';
import { Instance } from '@/components/stats/types';
import { ThemeData as ThemeProjectBreakdownData } from '@/components/stats/theme-project-breakdown-chart';
import { sectionMargins } from '@/theme';
import ThemeShelf from '@/components/themes/theme-shelf';
import { ThemeData } from '@/components/themes/theme-book';

import oscc1Fallback from '/public/images/hero/oscc1-fallback.png'
import xeniumFallback from '/public/images/hero/xenium-fallback.png'
import tonsilFallback from '/public/images/hero/tonsil-fallback.png'
import aboutImg1 from '/public/images/about-carousel/carousel-img-1.png'
import aboutImg2 from '/public/images/about-carousel/carousel-img-2.png'
import aboutImg3 from '/public/images/about-carousel/carousel-img-3.png'
import aboutImg4 from '/public/images/about-carousel/carousel-img-4.png'
import autoimmunityThemeImg from '/public/images/themes/Autoimmunity.png'
import cancerThemeImg from '/public/images/themes/Cancer.png'
import earlyLifeThemeImg from '/public/images/themes/Early Life.png'
import fibrosisThemeImg from '/public/images/themes/Fibrosis.png'
import healthyReferenceThemeImg from '/public/images/themes/Healthy Reference.png'
import infectionThemeImg from '/public/images/themes/Infection.png'
import inflammationThemeImg from '/public/images/themes/Inflammation.png'
import neurodegenerationThemeImg from '/public/images/themes/Neurodegeneration.png'
import womensHealthThemeImg from '/public/images/themes/Women\'s Health.png'


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
const daysInRange = 365 * 4
startDate.setDate(startDate.getDate() - daysInRange)
const dayInMs = 24 * 60 * 60 * 1000

const STATS = {
  bytes: Array(daysInRange).fill(null).map((_, i) => ({
    date: new Date(startDate.getTime() + i * dayInMs),
    value: getRandomArbitrary(1e12, 999e12),
  })),
  assays: Array(daysInRange).fill(null).map((_, i) => ({
    date: new Date(startDate.getTime() + i * dayInMs),
    value: getRandomArbitrary(1000, 400000),
  })),
  subjects: Array(daysInRange).fill(null).map((_, i) => ({
    date: new Date(startDate.getTime() + i * dayInMs),
    value: getRandomArbitrary(10, 999),
  })),
  files: Array(daysInRange).fill(null).map((_, i) => ({
    date: new Date(startDate.getTime() + i * dayInMs),
    value: getRandomArbitrary(1000, 999999),
  })),
  samples: Array(daysInRange).fill(null).map((_, i) => ({
    date: new Date(startDate.getTime() + i * dayInMs),
    value: getRandomArbitrary(100000, 999999),
  })),
  users: Array(daysInRange).fill(null).map((_, i) => ({
    date: new Date(startDate.getTime() + i * dayInMs),
    value: getRandomArbitrary(100, 1000),
  })),
}


const THEMES: ThemeData[] = [
  {
    name: 'Infection',
    color: '#89A7CE',
    project_count: 11,
    projects_link: '',
    description: faker.commerce.productDescription(),
    image: infectionThemeImg
  },
  {
    name: 'Autoimmunity',
    color: '#DDA373',
    project_count: 23,
    projects_link: '',
    description: faker.commerce.productDescription(),
    image: autoimmunityThemeImg
  },
  {
    name: 'Inflammation',
    color: '#D6D8A8',
    project_count: 8,
    projects_link: '',
    description: faker.commerce.productDescription(),
    image: inflammationThemeImg
  },
  {
    name: 'Fibrosis',
    color: '#A2A648',
    project_count: 2,
    projects_link: '',
    description: faker.commerce.productDescription(),
    image: fibrosisThemeImg
  },
  {
    name: 'Early Life'
    , color: '#7FA190'
    , project_count: 8,
    projects_link: '',
    description: faker.commerce.productDescription(),
    image: earlyLifeThemeImg
  },
  {
    name: 'Cancer',
    color: '#E4B8C7',
    project_count: 21,
    projects_link: '',
    description: faker.commerce.productDescription(),
    image: cancerThemeImg
  },
  {
    name: 'Neurodegeneration',
    color: '#556E66',
    project_count: 6,
    projects_link: '',
    description: faker.commerce.productDescription(),
    image: neurodegenerationThemeImg
  },
  {
    name: "Women's Health",
    color: '#E9C54E',
    project_count: 10,
    projects_link: '',
    description: faker.commerce.productDescription(),
    image: womensHealthThemeImg
  },
  {
    name: 'Healthy Reference',
    color: '#DFDED6',
    project_count: 4,
    projects_link: '',
    description: faker.commerce.productDescription(),
    image: healthyReferenceThemeImg
  },
]

THEMES.forEach(theme => {
  theme.projects_link = `#projects__${_.kebabCase(theme.name)}`
})

interface Project {

}

const PROJECTS: Project[] = []
Object.entries(THEMES).forEach(([k, v]) => {
  for (let i = 0; i < v.project_count; i++) {
    const project: Project = {
      name: faker.company.name(),
      theme: k,
    }

    PROJECTS.push(project)
  }
})

async function getData() {
  return {
    heroVideo: getRandomItem(VIDEOS),
    // TODO: replace with real data
    stats: STATS,
    // TODO: replace with real data
    aboutItems: ABOUT_ITEMS,
    themes: THEMES,
    projects: PROJECTS,
  }
}

export default async function Home() {
  const data = await getData()

  const carouselStats = {} as Record<keyof typeof STATS, number>
  for (const [k, v] of (Object.entries(STATS) as [keyof typeof STATS, Instance<number>[]][])) {
    carouselStats[k] = v[v.length - 1].value
  }

  const themeProjectBreakdown: ThemeProjectBreakdownData[] = []
  data.themes.forEach((theme) => {
    themeProjectBreakdown.push({
      name: theme.name,
      color: theme.color,
      project_count: theme.project_count,
    })
  })

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
            themeProjectBreakdown={themeProjectBreakdown}
          />
        </Box>

        <Box id='themes'>
          <ThemeShelf
            themeData={data.themes}
          />
        </Box>
      </Box>

    </React.Fragment>
  );
}
