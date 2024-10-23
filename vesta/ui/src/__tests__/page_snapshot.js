import { render } from '@testing-library/react'
import {VIDEOS, ABOUT_ITEMS} from '../app/page'

import { AppRouterCacheProvider } from '@mui/material-nextjs/v14-appRouter';
import { ThemeProvider } from '@mui/material/styles';
import theme from '@/theme';

import Hero from '@/components/hero/hero';


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

jest.mock('next/navigation', () => {
  return {
    __esModule: true,
    usePathname: () => ({
      pathname: '',
    }),
    useRouter: () => ({
      push: jest.fn(),
      replace: jest.fn(),
      prefetch: jest.fn(),
    }),
    useSearchParams: jest.fn(() => new URLSearchParams('param1=value1')),
    useServerInsertedHTML: jest.fn()
  }
})

const mockThemes = [
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
]

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

  it('ThemeShelf', () => {
      const { container } = render(
        <html lang="en">
          <body>
            <AppRouterCacheProvider>
              <ThemeProvider theme={theme}>
                <main>
                  <ThemeShelf
                    themeData={mockThemes}
                  />
                </main>
              </ThemeProvider>
            </AppRouterCacheProvider>
          </body>
        </html>)
      expect(container).toMatchSnapshot()
    })
})