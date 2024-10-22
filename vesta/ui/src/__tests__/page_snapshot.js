import { render, screen } from '@testing-library/react'
import {VIDEOS} from '../app/page'

import { AppRouterCacheProvider } from '@mui/material-nextjs/v14-appRouter';
import { ThemeProvider } from '@mui/material/styles';
import theme from '@/theme';

import Hero from '@/components/hero/hero';

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
})