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
import { ThemeProjectBreakdownData } from '@/components/stats/theme-project-breakdown-chart';
import { sectionMargins } from '@/theme';
import ThemeShelf from '@/components/themes/theme-shelf';
import { ThemeData } from '@/components/themes/models';
import ProjectListings from '@/components/project-listings/project-listings';
import { Project, ProjectStatus, ProjectType } from '@/components/project-listings/models';

import oscc1Fallback from '/public/images/hero/oscc1-fallback.png'
import xeniumFallback from '/public/images/hero/xenium-fallback.png'
import tonsilFallback from '/public/images/hero/tonsil-fallback.png'
import aboutImg1 from '/public/images/about-carousel/carousel-img-1.png'
import aboutImg2 from '/public/images/about-carousel/carousel-img-2.png'
import aboutImg3 from '/public/images/about-carousel/carousel-img-3.png'
import aboutImg4 from '/public/images/about-carousel/carousel-img-4.png'

import autoimmunityThemeImg from '/public/images/themes/covers/autoimmunity.png'
import cancerThemeImg from '/public/images/themes/covers/cancer.png'
import earlyLifeThemeImg from '/public/images/themes/covers/early-life.png'
import fibrosisThemeImg from '/public/images/themes/covers/fibrosis.png'
import healthyReferenceThemeImg from '/public/images/themes/covers/healthy-reference.png'
import infectionThemeImg from '/public/images/themes/covers/infection.png'
import inflammationThemeImg from '/public/images/themes/covers/inflammation.png'
import neurodegenerationThemeImg from '/public/images/themes/covers/neurodegeneration.png'
import womensHealthThemeImg from '/public/images/themes/covers/womens-health.png'

import autoimmunityThemeIcon from '/public/images/themes/icons/lg/svg/autoimmunity.svg'
import cancerThemeIcon from '/public/images/themes/icons/lg/svg/cancer.svg'
import earlyLifeThemeIcon from '/public/images/themes/icons/lg/svg/early-life.svg'
import fibrosisThemeIcon from '/public/images/themes/icons/lg/svg/fibrosis.svg'
import healthyReferenceThemeIcon from '/public/images/themes/icons/lg/svg/healthy-reference.svg'
import infectionThemeIcon from '/public/images/themes/icons/lg/svg/infection.svg'
import inflammationThemeIcon from '/public/images/themes/icons/lg/svg/inflammation.svg'
import neurodegenerationThemeIcon from '/public/images/themes/icons/lg/svg/neurodegeneration.svg'
import womensHealthThemeIcon from '/public/images/themes/icons/lg/svg/womens-health.svg'

import autoimmunityThemeFiltered from '/public/images/themes/components/autoimmunity-filteredimg.png'
import cancerThemeFiltered from '/public/images/themes/components/cancer-filteredimg.png'
import earlyLifeThemeFiltered from '/public/images/themes/components/early-life-filteredimg.png'
import fibrosisThemeFiltered from '/public/images/themes/components/fibrosis-filteredimg.png'
import healthyReferenceThemeFiltered from '/public/images/themes/components/healthy-reference-filteredimg.png'
import infectionThemeFiltered from '/public/images/themes/components/infection-filteredimg.png'
import inflammationThemeFiltered from '/public/images/themes/components/inflammation-filteredimg.png'
import neurodegenerationThemeFiltered from '/public/images/themes/components/neurodegeneration-filteredimg.png'
import womensHealthThemeFiltered from '/public/images/themes/components/womens-health-filteredimg.png'

import autoimmunityThemeProjectBg from '/public/images/themes/components/autoimmunity-projbg.png'
import cancerThemeProjectBg from '/public/images/themes/components/cancer-projbg.png'
import earlyLifeThemeProjectBg from '/public/images/themes/components/early-life-projbg.png'
import fibrosisThemeProjectBg from '/public/images/themes/components/fibrosis-projbg.png'
import healthyReferenceThemeProjectBg from '/public/images/themes/components/healthy-reference-projbg.png'
import infectionThemeProjectBg from '/public/images/themes/components/infection-projbg.png'
import inflammationThemeProjectBg from '/public/images/themes/components/inflammation-projbg.png'
import neurodegenerationThemeProjectBg from '/public/images/themes/components/neurodegeneration-projbg.png'
import womensHealthThemeProjectBg from '/public/images/themes/components/womens-health-projbg.png'


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
debugger
const startDate = new Date(Date.now())
const daysInRange = 365 * 4
startDate.setDate(startDate.getDate() - daysInRange)
const dayInMs = 24 * 60 * 60 * 1000

const STATS = {
  bytes: Array(daysInRange).fill(null).map((_, i) => ({
    date: new Date(startDate.getTime() + i * dayInMs),
    value: faker.helpers.rangeToNumber({ min: 1e12, max: 999e12, }),
  })),
  assays: Array(daysInRange).fill(null).map((_, i) => ({
    date: new Date(startDate.getTime() + i * dayInMs),
    value: faker.helpers.rangeToNumber({ min: 1000, max: 400000, }),
  })),
  subjects: Array(daysInRange).fill(null).map((_, i) => ({
    date: new Date(startDate.getTime() + i * dayInMs),
    value: faker.helpers.rangeToNumber({ min: 10, max: 999, }),
  })),
  files: Array(daysInRange).fill(null).map((_, i) => ({
    date: new Date(startDate.getTime() + i * dayInMs),
    value: faker.helpers.rangeToNumber({ min: 1000, max: 999999, }),
  })),
  samples: Array(daysInRange).fill(null).map((_, i) => ({
    date: new Date(startDate.getTime() + i * dayInMs),
    value: faker.helpers.rangeToNumber({ min: 100000, max: 999999, }),
  })),
  users: Array(daysInRange).fill(null).map((_, i) => ({
    date: new Date(startDate.getTime() + i * dayInMs),
    value: faker.helpers.rangeToNumber({ min: 100, max: 1000, }),
  })),
}


const THEMES: ThemeData[] = [
  {
    name: 'Infection',
    color: '#89A7CE',
    textColor: 'dark',
    projectCount: 11,
    projectsLink: '',
    description: faker.commerce.productDescription(),
    imageComponents: {
      filtered: infectionThemeFiltered,
      projectBackground: infectionThemeProjectBg,
    },
    coverImage: infectionThemeImg,
    icon: infectionThemeIcon,
  },
  {
    name: 'Autoimmunity',
    color: '#DDA373',
    textColor: 'dark',
    projectCount: 23,
    projectsLink: '',
    description: faker.commerce.productDescription(),
    imageComponents: {
      filtered: autoimmunityThemeFiltered,
      projectBackground: autoimmunityThemeProjectBg,
    },
    coverImage: autoimmunityThemeImg,
    icon: autoimmunityThemeIcon,
  },
  {
    name: 'Inflammation',
    color: '#D6D8A8',
    textColor: 'dark',
    projectCount: 8,
    projectsLink: '',
    description: faker.commerce.productDescription(),
    imageComponents: {
      filtered: inflammationThemeFiltered,
      projectBackground: inflammationThemeProjectBg,
    },
    coverImage: inflammationThemeImg,
    icon: inflammationThemeIcon,
  },
  {
    name: 'Fibrosis',
    color: '#A2A648',
    textColor: 'dark',
    projectCount: 2,
    projectsLink: '',
    description: faker.commerce.productDescription(),
    imageComponents: {
      filtered: fibrosisThemeFiltered,
      projectBackground: fibrosisThemeProjectBg,
    },
    coverImage: fibrosisThemeImg,
    icon: fibrosisThemeIcon,
  },
  {
    name: 'Early Life',
    color: '#7FA190',
    textColor: 'dark',
    projectCount: 8,
    projectsLink: '',
    description: faker.commerce.productDescription(),
    imageComponents: {
      filtered: earlyLifeThemeFiltered,
      projectBackground: earlyLifeThemeProjectBg,
    },
    coverImage: earlyLifeThemeImg,
    icon: earlyLifeThemeIcon,
  },
  {
    name: 'Cancer',
    color: '#E4B8C7',
    textColor: 'dark',
    projectCount: 21,
    projectsLink: '',
    description: faker.commerce.productDescription(),
    imageComponents: {
      filtered: cancerThemeFiltered,
      projectBackground: cancerThemeProjectBg,
    },
    coverImage: cancerThemeImg,
    icon: cancerThemeIcon,
  },
  {
    name: 'Neurodegeneration',
    color: '#556E66',
    textColor: 'light',
    projectCount: 6,
    projectsLink: '',
    description: faker.commerce.productDescription(),
    imageComponents: {
      filtered: neurodegenerationThemeFiltered,
      projectBackground: neurodegenerationThemeProjectBg,
    },
    coverImage: neurodegenerationThemeImg,
    icon: neurodegenerationThemeIcon,
  },
  {
    name: "Women's Health",
    color: '#E9C54E',
    textColor: 'dark',
    projectCount: 10,
    projectsLink: '',
    description: faker.commerce.productDescription(),
    imageComponents: {
      filtered: womensHealthThemeFiltered,
      projectBackground: womensHealthThemeProjectBg,
    },
    coverImage: womensHealthThemeImg,
    icon: womensHealthThemeIcon,
  },
  {
    name: 'Healthy Reference',
    color: '#DFDED6',
    textColor: 'dark',
    projectCount: 4,
    projectsLink: '',
    description: faker.commerce.productDescription(),
    imageComponents: {
      filtered: healthyReferenceThemeFiltered,
      projectBackground: healthyReferenceThemeProjectBg,
    },
    coverImage: healthyReferenceThemeImg,
    icon: healthyReferenceThemeIcon,
  },
]

THEMES.forEach(theme => {
  theme.projectsLink = `#projects__${_.kebabCase(theme.name)}`
})

function createRandomArray<E extends any>(min: number, max: number, generator: () => E): Array<E> {
  const result = [] as Array<E>
  for (let i = 0; i < faker.helpers.rangeToNumber({ min: min, max: max, }); i++) {
    result.push(generator())
  }
  return result
}

const PROJECTS: Project[] = []
THEMES.forEach(theme => {
  for (let i = 0; i < theme.projectCount; i++) {
    const name = faker.company.name().replace(' -', '')

    const project: Project = {
      name: name.toLowerCase().slice(0, 4),
      fullName: name,
      heading: Math.round(Math.random()) === 1 ? faker.lorem.sentence() : undefined,
      description: faker.lorem.paragraphs({ min: 1, max: 3 }),
      fundingSource: faker.commerce.department(),
      principalInvestigators: createRandomArray(1, 5, () => faker.person.fullName()),
      status: faker.helpers.arrayElement(Object.values(ProjectStatus)),
      type: faker.helpers.arrayElement(Object.values(ProjectType)),
      species: faker.lorem.word(),
      startDate: faker.date.between({ from: '2018-01-01T00:00:00.000Z', to: '2024-01-01T00:00:00.000Z' }),
      dataCollectionComplete: Math.round(Math.random()) === 1,
      userCount: faker.helpers.rangeToNumber({ min: 5, max: 100, }),
      theme: theme,
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
      projectCount: theme.projectCount,
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

        <Box id="projects">
          <ProjectListings
            projectData={_.sortBy(data.projects, (p) => p.fullName)}
          />
        </Box>
      </Box>

    </React.Fragment>
  );
}
