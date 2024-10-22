import { render, screen } from '@testing-library/react'
import Home, {VIDEOS, ABOUT_ITEMS} from '../app/page'
import { faker } from '@faker-js/faker'

import { AppRouterCacheProvider } from '@mui/material-nextjs/v14-appRouter';
import { ThemeProvider } from '@mui/material/styles';
import theme from '@/theme';

import Hero from '@/components/hero/hero';

import ProjectListings from '@/components/project-listings/project-listings';
import _ from 'lodash'

import AboutCarousel from '@/components/about/about-carousel';

import ThemeShelf from '@/components/themes/theme-shelf';

// Mock Swiper components (https://stackoverflow.com/a/70685443)
jest.mock('swiper/react', () => ({
  Swiper: ({ children }) => (
    <div data-testid="swiper-testid">{children}</div>
  ),
  SwiperSlide: ({ children }) => (
    <div data-testid="swiper-slide-testid">{children}</div>
  ),
  SwiperClass: ({ children }) => (
    <div data-testid="swiper-class-testid">{children}</div>
  ),
}))
jest.mock('swiper/modules', () => ({
  EffectFade: (props) => null,
  Autoplay: (props) => null,
  A11y: (props) => null,
}))

// // Mock the `getData` function element
// jest.mock('../app/page', () => {
//   const originalModule = jest.requireActual('../app/page');
//   const mockData = {
//     heroVideos: originalModule.VIDEOS,
//     stats: {
//       global: {
//         sampleCount: 1000,
//         assayCount: 100,
//         userCount: 50,
//       }
//     },
//     aboutItems: originalModule.ABOUT_ITEMS,
//     themes: [
//       {
//         name: 'Infection',
//         color: '#89A7CE',
//         altColor: '#E8F2FF',
//         textColor: 'light',
//         projectCount: 1,
//         projectsLink: '',
//         description: 'Infection is bad',
//         imageComponents: {
//           filtered: jest.requireActual('/public/images/themes/components/infection-filteredimg.png'),
//           projectBackground: jest.requireActual('/public/images/themes/components/infection-projbg.png'),
//         },
//         coverImage: jest.requireActual('/public/images/themes/covers/infection.png'),
//         icon: jest.requireActual('/public/images/themes/icons/infection.svg'),
//       },
//       {
//         name: 'Cancer',
//         color: '#E4B8C7',
//         altColor: '#B53B38',
//         textColor: 'dark',
//         projectCount: 2,
//         projectsLink: '',
//         description: 'cancer lorem ipsum',
//         imageComponents: {
//           filtered: jest.requireActual('/public/images/themes/components/cancer-filteredimg.png'),
//           projectBackground: jest.requireActual('/public/images/themes/components/cancer-projbg.png'),
//         },
//         coverImage: jest.requireActual('/public/images/themes/covers/cancer.png'),
//         icon: jest.requireActual('/public/images/themes/icons/cancer.svg'),
//       }
//     ],
//     projects: [
//       {
//         name: "labours",
//         fullName: "The Labours Project",
//         description: "lorem ipsum",
//         fundingSource: "lorem ipsum",
//         principalInvestigators: [
//           {
//             name: "John Doe",
//             email: "j.d@email.com",
//             title: "Title",
//             imageUrl: undefined,
//             profileUrl: undefined,
//             color: '#E4B8C7',
//             altColor: 'utilityHighlight.main',
//           }
//         ],
//         status: "Team",
//         type: "CoProject",
//         dataTypes: ['clinical', 'scRNAseq'],
//         sampleCount: 27,
//         assayCount: 84,
//         hasClinicalData: 'Yes',
//         species: "Human",
//         startDate: "2020",
//         dataCollectionComplete: true,
//         userCount: 4,
//         theme: "Cancer",
//         href: (new URL('/labours', 'https://example.com')).href
//       },
//       {
//         name: "victims",
//         fullName: "The Victims Project",
//         description: "lorem ipsum",
//         fundingSource: "lorem ipsum",
//         principalInvestigators: [
//           {
//             name: "Jane Doe",
//             email: "j.d@email.com",
//             title: "Title",
//             imageUrl: undefined,
//             profileUrl: undefined,
//             color: '#E4B8C7',
//             altColor: 'utilityHighlight.main',
//           },
//           {
//             name: "John Doe",
//             email: "j.d@email.com",
//             title: "Title",
//             imageUrl: undefined,
//             profileUrl: undefined,
//             color: '#E4B8C7',
//             altColor: 'utilityHighlight.main',
//           }
//         ],
//         status: "Team",
//         type: "CoProject",
//         dataTypes: ['scRNAseq', 'bulkRNAseq', 'treatment'],
//         sampleCount: 48,
//         assayCount: 312,
//         hasClinicalData: 'Yes',
//         species: "Human",
//         startDate: "2022",
//         dataCollectionComplete: true,
//         userCount: 8,
//         theme: "Infection",
//         href: (new URL('/victims', 'https://example.com')).href
//       }
//     ],
//     accessUrl: 'https://example.com'
//   };

//   return {
//     __esModule: true,
//     ...originalModule,
//     getData: jest.fn(async () => await mockData),
//     Home
//   };
// });

const mockData = {
  heroVideos: VIDEOS,
  stats: {
    global: {
      1698002688000: {
        sampleCount: 1000,
        assayCount: 100,
        userCount: 50,
      },
      1698002687000: {
        sampleCount: 900,
        assayCount: 90,
        userCount: 45,
      }
    }
  },
  aboutItems: ABOUT_ITEMS,
  themes: [
    {
      name: 'Infection',
      color: '#89A7CE',
      altColor: '#E8F2FF',
      textColor: 'light',
      projectCount: 1,
      projectsLink: '',
      description: 'Infection is bad',
      imageComponents: {
        filtered: jest.requireActual('/public/images/themes/components/infection-filteredimg.png'),
        projectBackground: jest.requireActual('/public/images/themes/components/infection-projbg.png'),
      },
      coverImage: jest.requireActual('/public/images/themes/covers/infection.png'),
      icon: jest.requireActual('/public/images/themes/icons/infection.svg'),
    },
    {
      name: 'Cancer',
      color: '#E4B8C7',
      altColor: '#B53B38',
      textColor: 'dark',
      projectCount: 2,
      projectsLink: '',
      description: 'cancer lorem ipsum',
      imageComponents: {
        filtered: jest.requireActual('/public/images/themes/components/cancer-filteredimg.png'),
        projectBackground: jest.requireActual('/public/images/themes/components/cancer-projbg.png'),
      },
      coverImage: jest.requireActual('/public/images/themes/covers/cancer.png'),
      icon: jest.requireActual('/public/images/themes/icons/cancer.svg'),
    }
  ],
  projects: [
    {
      name: "labours",
      fullName: "The Labours Project",
      description: "lorem ipsum",
      fundingSource: "lorem ipsum",
      principalInvestigators: [
        {
          name: "John Doe",
          email: "j.d@email.com",
          title: "Title",
          imageUrl: undefined,
          profileUrl: undefined,
          color: '#E4B8C7',
          altColor: 'utilityHighlight.main',
        }
      ],
      status: "Team",
      type: "CoProject",
      dataTypes: ['clinical', 'scRNAseq'],
      sampleCount: 27,
      assayCount: 84,
      hasClinicalData: 'Yes',
      species: "Human",
      startDate: "2020",
      dataCollectionComplete: true,
      userCount: 4,
      theme: "Cancer",
      href: (new URL('/labours', 'https://example.com')).href
    },
    {
      name: "victims",
      fullName: "The Victims Project",
      description: "lorem ipsum",
      fundingSource: "lorem ipsum",
      principalInvestigators: [
        {
          name: "Jane Doe",
          email: "j.d@email.com",
          title: "Title",
          imageUrl: undefined,
          profileUrl: undefined,
          color: '#E4B8C7',
          altColor: 'utilityHighlight.main',
        },
        {
          name: "John Doe",
          email: "j.d@email.com",
          title: "Title",
          imageUrl: undefined,
          profileUrl: undefined,
          color: '#E4B8C7',
          altColor: 'utilityHighlight.main',
        }
      ],
      status: "Team",
      type: "CoProject",
      dataTypes: ['scRNAseq', 'bulkRNAseq', 'treatment'],
      sampleCount: 48,
      assayCount: 312,
      hasClinicalData: 'Yes',
      species: "Human",
      startDate: "2022",
      dataCollectionComplete: true,
      userCount: 8,
      theme: "Infection",
      href: (new URL('/victims', 'https://example.com')).href
    }
  ],
  accessUrl: 'https://example.com'
};

describe('Homepage', () => {
  beforeEach(() => {
    // Mocking faker.helpers.rangeToNumber to always return 0
    jest.spyOn(faker.helpers, 'rangeToNumber').mockReturnValue(0);
    jest.spyOn(require('../app/page'), "getData").mockReturnValue(mockData);
    global.fetch = jest.fn(() => Promise.resolve({
      json: () => Promise.resolve({ data: 'some data' })
    }));
  });

  afterEach(() => {
    // Restore the original implementation after each test
    jest.restoreAllMocks();
  });


  xit('renders homepage unchanged', async () => {
    const { container } = render(<Home/>)
    await screen.findByText('About the Library');
    expect(container).toMatchSnapshot()
  })
})

describe('Renderable homepage elements', () => {
  beforeEach(()=>{
    jest.spyOn(require('next/router'), 'useRouter')
  })

  afterEach(() => {
    // Restore the original implementation after each test
    jest.restoreAllMocks();
  });

  it('Hero', () => {
    const { container } = render(
      <html lang="en">
        <body>
          <AppRouterCacheProvider>
            <ThemeProvider theme={theme}>
              <main>
                <Hero
                  videos={VIDEOS}
                  initVideoIdx={0}
                  stats={{
                    bytes: 1000005,
                    assays: 100,
                    subjects: 148,
                    files: 597,
                    samples: 1000,
                    users: 50,
                  }}
                  scrollTargetId='about'
                  accessUrl={'https://example.com'}
                />
              </main>
            </ThemeProvider>
          </AppRouterCacheProvider>
        </body>
      </html>)
    expect(container).toMatchSnapshot()
  })

  xit('ThemeShelf', () => {
    const { container } = render(
      <html lang="en">
        <body>
          <AppRouterCacheProvider>
            <ThemeProvider theme={theme}>
              <main>
                <ThemeShelf
                  themeData={mockData.themes}
                />
              </main>
            </ThemeProvider>
          </AppRouterCacheProvider>
        </body>
      </html>)
    expect(container).toMatchSnapshot()
  })

  xit('ProjectListings', () => {
    const { container } = render(
      <html lang="en">
        <body>
          <AppRouterCacheProvider>
            <ThemeProvider theme={theme}>
              <main>
                <ProjectListings
                  projectData={_.sortBy(mockData.projects, (p) => p.fullName)}
                />
              </main>
            </ThemeProvider>
          </AppRouterCacheProvider>
        </body>
      </html>)
    expect(container).toMatchSnapshot()
  })

  xit('AboutCarousel', () => {
    const { container } = render(
      <html lang="en">
        <head>
        </head>
        <body>
          <AppRouterCacheProvider>
            <ThemeProvider theme={theme}>
              {/* <CssBaseline enableColorScheme={true} />
              <MainNav
                accessUrl={'https://example.com'}
              /> */}
              <main>
                <AboutCarousel
                  items={ABOUT_ITEMS}
                />
              </main>
            </ThemeProvider>
          </AppRouterCacheProvider>
        </body>
      </html>
    )
    expect(container).toMatchSnapshot()

    const element = getByText('About the Library');
    expect(element).toBeInTheDocument();
  })
})

describe('Renderable homepage elements', () => {
  
})