import * as React from 'react';
import { Metadata } from 'next';
import Box from '@mui/system/Box'
import { faker } from '@faker-js/faker'
import _ from 'lodash'

import Hero from '@/components/hero/hero';
import AboutCarousel from '@/components/about/about-carousel';
import LibraryStats from '@/components/stats/library-stats';
import { Instance } from '@/components/stats/types';
import { ThemeProjectBreakdownData } from '@/components/stats/theme-project-breakdown-chart';
import { sectionMargins } from '@/theme';
import ThemeShelf from '@/components/themes/theme-shelf';
import { ThemeData } from '@/components/themes/models';
import ProjectListings from '@/components/project-listings/project-listings';
import { Project, ProjectStatus, ProjectType, ProjectsSearchParamsState, PROJECTS_SEARCH_PARAMS_KEY } from '@/components/project-listings/models';
import { toSearchParamsString } from '@/lib/utils/uri';
import { VestaApiClient } from '@/lib/clients/vesta-api/client';
import { defaultDict } from '@/lib/utils/object';

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

import autoimmunityThemeIcon from '/public/images/themes/icons/autoimmunity.svg'
import cancerThemeIcon from '/public/images/themes/icons/cancer.svg'
import earlyLifeThemeIcon from '/public/images/themes/icons/early-life.svg'
import fibrosisThemeIcon from '/public/images/themes/icons/fibrosis.svg'
import healthyReferenceThemeIcon from '/public/images/themes/icons/healthy-reference.svg'
import infectionThemeIcon from '/public/images/themes/icons/infection.svg'
import inflammationThemeIcon from '/public/images/themes/icons/inflammation.svg'
import neurodegenerationThemeIcon from '/public/images/themes/icons/neurodegeneration.svg'
import womensHealthThemeIcon from '/public/images/themes/icons/womens-health.svg'

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


// Cache configuration
// https://nextjs.org/docs/app/api-reference/file-conventions/route-segment-config#revalidate
export const revalidate = 86400  // 1 day

export const metadata: Metadata = {
  title: 'UCSF Data Library',
}

const VIDEOS = [
  { videoSrc: '/videos/hero/oscc1-clipped.mp4', imageSrc: oscc1Fallback },
  { videoSrc: '/videos/hero/xenium-clipped.mp4', imageSrc: xeniumFallback },
  { videoSrc: '/videos/hero/tonsil-clipped.mp4', imageSrc: tonsilFallback },
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
    title: 'Contribute Data',
    header: 'Contributing to the library',
    body: 'Designed as a resource for open collaboration, the UCSF Data library aims to be a catalyst for new research discoveries. We encourage everyone, both UCSF members and external collaborators, to share their data on our platform.',
    link: {
      header: 'Working on something great?',
      blurb: 'We’re always looking for new data and projects to make more accessible.',
      href: '/#contribute',
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
      href: process.env.TIMUR_URL,
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
      href: process.env.TIMUR_URL,
      label: 'Get Access',
    },
    image: {
      src: aboutImg4,
      alt: 'A few Data Library library cards',
    },
  },
]

const THEMES: ThemeData[] = [
  {
    name: 'Infection',
    color: '#89A7CE',
    altColor: '#E8F2FF',
    textColor: 'light',
    projectCount: 0,
    projectsLink: '',
    description: "Microbial infections, including viral and bacterial pathogenesis",
    imageComponents: {
      filtered: infectionThemeFiltered,
      projectBackground: infectionThemeProjectBg,
    },
    coverImage: infectionThemeImg,
    icon: infectionThemeIcon,
  },
  {
    name: 'Autoimmunity',
    color: '#D36F49',
    altColor: '#556E66',
    textColor: 'light',
    projectCount: 0,
    projectsLink: '',
    description: "Autoimmune and rheumatic conditions",
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
    altColor: '#F3F2E3',
    textColor: 'dark',
    projectCount: 0,
    projectsLink: '',
    description: "Conditions that confer a predominantly inflammatory phenotype",
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
    altColor: '#13283F',
    textColor: 'light',
    projectCount: 0,
    projectsLink: '',
    description: "Conditions that result in significant tissue fibrosis",
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
    altColor: '#3A4B48',
    textColor: 'light',
    projectCount: 0,
    projectsLink: '',
    description: "Immune and systems development in early life",
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
    altColor: '#B53B38',
    textColor: 'dark',
    projectCount: 0,
    projectsLink: '',
    description: "Human cancers and animal models, including affected tissue and peripheral sampling",
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
    altColor: '#2A5A8D',
    textColor: 'light',
    projectCount: 0,
    projectsLink: '',
    description: "Conditions affecting the brain and nervous system",
    imageComponents: {
      filtered: neurodegenerationThemeFiltered,
      projectBackground: neurodegenerationThemeProjectBg,
    },
    coverImage: neurodegenerationThemeImg,
    icon: neurodegenerationThemeIcon,
  },
  {
    name: "Womens Health",
    color: '#E9C54E',
    altColor: '#B53B38',
    textColor: 'dark',
    projectCount: 0,
    projectsLink: '',
    description: "Conditions predominantly affecting women such as those affecting reproductive organs",
    imageComponents: {
      filtered: womensHealthThemeFiltered,
      projectBackground: womensHealthThemeProjectBg,
    },
    coverImage: womensHealthThemeImg,
    icon: womensHealthThemeIcon,
  },
  {
    name: 'Healthy Reference',
    color: '#DDDAD0',
    altColor: '#F9F8F6',
    textColor: 'dark',
    projectCount: 0,
    projectsLink: '',
    description: "Studies of healthy individuals for use as reference in pathologic settings",
    imageComponents: {
      filtered: healthyReferenceThemeFiltered,
      projectBackground: healthyReferenceThemeProjectBg,
    },
    coverImage: healthyReferenceThemeImg,
    icon: healthyReferenceThemeIcon,
  },
]

// Manage projects link
THEMES.forEach(theme => {
  const projectsSearchParamsState: ProjectsSearchParamsState = {
    filters: { theme: [theme.name] }
  }
  const search = toSearchParamsString({ [PROJECTS_SEARCH_PARAMS_KEY]: projectsSearchParamsState })
  theme.projectsLink = '/?' + search + '#projects'
})


async function getData() {
  const apiClient = new VestaApiClient()

  const [
    stats,
    apiProjects,
  ] = await Promise.all([
    apiClient.fetchStats(),
    apiClient.fetchProjects(),
  ])

  const themeProjectsCount = defaultDict<string, number>(_ => 0)

  const projects: (Project | undefined)[] = apiProjects.map(proj => {
    const status = Object.values(ProjectStatus).find(k => k.toUpperCase() === proj.status.toUpperCase())
    const type = Object.values(ProjectType).find(k => k.toUpperCase() === proj.type.toUpperCase())
    const theme = THEMES.find(theme => theme.name.toUpperCase() === proj.theme.toUpperCase())
    if (!status || !type || !theme) {
      return
    }

    themeProjectsCount[theme.name.toUpperCase()] += 1

    const dataTypes = proj.data_types.filter(dt => dt.toUpperCase() !== 'project'.toUpperCase())

    return {
      name: proj.name,
      fullName: _.words(proj.full_name).join(' '),
      description: proj.description,
      fundingSource: proj.funding_source,
      principalInvestigators: proj.principal_investigators.map(pi => {
        const theme = faker.helpers.arrayElement(THEMES)

        return {
          name: pi.name,
          email: pi.email,
          title: pi.title,
          imageUrl: pi.image_url,
          profileUrl: pi.profile_url,
          color: theme.color,
          altColor: theme.textColor === 'light' ? 'utilityHighlight.main' : 'ground.grade10',
        }
      }),
      status,
      type,
      dataTypes,
      sampleCount: proj.name in stats.byProjectName && stats.byProjectName[proj.name].samples.length > 0 ? stats.byProjectName[proj.name].samples.at(-1)?.value || 0 : 0,
      assayCount: proj.name in stats.byProjectName && stats.byProjectName[proj.name].assays.length > 0 ? stats.byProjectName[proj.name].assays.at(-1)?.value || 0 : 0,
      hasClinicalData: (proj.name in stats.byProjectName && stats.byProjectName[proj.name].assays.length > 0 ? stats.byProjectName[proj.name].assays.at(-1)?.value || 0 : 0) > 0 ? 'Yes' : 'No',
      species: proj.species,
      startDate: proj.start_date,
      dataCollectionComplete: proj.data_collection_complete,
      userCount: proj.name in stats.byProjectName && stats.byProjectName[proj.name].users.length > 0 ? stats.byProjectName[proj.name].users.at(-1)?.value || 0 : 0,
      theme,
      href: (new URL(`/${proj.name}`, process.env.TIMUR_URL)).href
    }
  })

  const themes = THEMES.map(theme => ({
    ...theme,
    projectCount: themeProjectsCount[theme.name.toUpperCase()],
  }))

  return {
    heroVideos: VIDEOS,
    stats: stats.global,
    aboutItems: ABOUT_ITEMS,
    themes,
    projects: (projects.filter(p => p) as Project[]),
    accessUrl: process.env.TIMUR_URL
  }
}

export default async function Home() {
  const data = await getData()

  const carouselStats = {} as Record<keyof typeof data.stats, number>
  for (const [k, v] of (Object.entries(data.stats) as [keyof typeof data.stats, Instance<number>[]][])) {
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
        videos={data.heroVideos}
        initVideoIdx={faker.helpers.rangeToNumber({ min: 0, max: data.heroVideos.length - 1 })}
        stats={carouselStats}
        scrollTargetId='about'
        accessUrl={data.accessUrl}
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

        <Box id='projects'>
          <ProjectListings
            projectData={_.sortBy(data.projects, (p) => p.fullName)}
          />
        </Box>
      </Box>

    </React.Fragment>
  );
}
