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
import { getData } from '@/lib/clients/vesta-api/request';

// Cache configuration
// https://nextjs.org/docs/app/api-reference/file-conventions/route-segment-config#revalidate
export const revalidate = 86400  // 1 day

export const metadata: Metadata = {
  title: 'UCSF Data Library',
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
      </Box>

    </React.Fragment>
  );
}
