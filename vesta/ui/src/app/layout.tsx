import * as React from 'react';
import { AppRouterCacheProvider } from '@mui/material-nextjs/v14-appRouter';
import { ThemeProvider } from '@mui/material/styles';
import CssBaseline from '@mui/material/CssBaseline';
import 'swiper/css';
import 'swiper/css/bundle';
import { faker } from '@faker-js/faker';

import theme from '@/theme';
import MainNav from '@/components/nav/main-nav';
import Footer from '@/components/footer/footer';
import { User } from '@/components/user/models';
import { UserContextProvider } from '@/components/user/context';

import sample1 from '/public/images/library-card/image-samples/1.png'
import sample2 from '/public/images/library-card/image-samples/2.png'
import sample3 from '/public/images/library-card/image-samples/3.png'
import sample4 from '/public/images/library-card/image-samples/4.png'
import sample5 from '/public/images/library-card/image-samples/5.png'
import sample6 from '/public/images/library-card/image-samples/6.png'
import sample7 from '/public/images/library-card/image-samples/7.png'
import sample8 from '/public/images/library-card/image-samples/8.png'
import sample9 from '/public/images/library-card/image-samples/9.png'
import sample10 from '/public/images/library-card/image-samples/10.png'
import sample11 from '/public/images/library-card/image-samples/11.png'
import sample12 from '/public/images/library-card/image-samples/12.png'
import sample13 from '/public/images/library-card/image-samples/13.png'
import sample14 from '/public/images/library-card/image-samples/14.png'
import sample15 from '/public/images/library-card/image-samples/15.png'
import sample16 from '/public/images/library-card/image-samples/16.png'
import sample17 from '/public/images/library-card/image-samples/17.png'
import sample18 from '/public/images/library-card/image-samples/18.png'
import sample19 from '/public/images/library-card/image-samples/19.png'
import sample20 from '/public/images/library-card/image-samples/20.png'
import sample21 from '/public/images/library-card/image-samples/21.png'
import sample22 from '/public/images/library-card/image-samples/22.png'
import sample23 from '/public/images/library-card/image-samples/23.png'
import sample24 from '/public/images/library-card/image-samples/24.png'
import sample25 from '/public/images/library-card/image-samples/25.png'
import sample26 from '/public/images/library-card/image-samples/26.png'
import sample27 from '/public/images/library-card/image-samples/27.png'
import sample28 from '/public/images/library-card/image-samples/28.png'
import sample29 from '/public/images/library-card/image-samples/29.png'
import sample30 from '/public/images/library-card/image-samples/30.png'
import sample31 from '/public/images/library-card/image-samples/31.png'
import sample32 from '/public/images/library-card/image-samples/32.png'
import sample33 from '/public/images/library-card/image-samples/33.png'
import sample34 from '/public/images/library-card/image-samples/34.png'
import sample35 from '/public/images/library-card/image-samples/35.png'
import sample36 from '/public/images/library-card/image-samples/36.png'

const IMAGES = [
  sample1,
  sample2,
  sample3,
  sample4,
  sample5,
  sample6,
  sample7,
  sample8,
  sample9,
  sample10,
  sample11,
  sample12,
  sample13,
  sample14,
  sample15,
  sample16,
  sample17,
  sample18,
  sample19,
  sample20,
  sample21,
  sample22,
  sample23,
  sample24,
  sample25,
  sample26,
  sample27,
  sample28,
  sample29,
  sample30,
  sample31,
  sample32,
  sample33,
  sample34,
  sample35,
  sample36,
]


const USER: User = {
  name: faker.person.fullName(),
  email: faker.internet.email(),
  title: faker.person.jobTitle(),
  role: faker.person.jobType(),
  imageUrl: faker.helpers.arrayElement(IMAGES),
  avatarUrl: faker.helpers.maybe(() => faker.image.avatar()),
  color: faker.color.rgb(),
  joinDate: faker.date.past({ years: 10 }),
  contributions: faker.number.int({ min: 0, max: 100 }),
  projectMemberships: faker.number.int({ min: 0, max: 100 }),
}

async function getData() {
  return {
    user: USER,
  }
}


export default async function RootLayout(props: { children: React.ReactNode }) {
  const data = await getData()

  return (
    <html lang="en">
      <head>
      </head>
      <body>
        <AppRouterCacheProvider>
          <ThemeProvider theme={theme}>
            <CssBaseline enableColorScheme={true} />

            <UserContextProvider
              value={data.user}
            >

              <MainNav />
              <main>{props.children}</main>
              <Footer />

            </UserContextProvider>

          </ThemeProvider>
        </AppRouterCacheProvider>
      </body>
    </html>
  );
}
