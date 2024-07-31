import * as React from 'react';
import { AppRouterCacheProvider } from '@mui/material-nextjs/v14-appRouter';
import { ThemeProvider } from '@mui/material/styles';
import CssBaseline from '@mui/material/CssBaseline';
import 'swiper/css';
import 'swiper/css/bundle';
import { faker } from '@faker-js/faker';

import theme from '@/theme';
import MainNav from '@/components/nav/main-nav';
import Footer from '@/components/footer';
import { User } from '@/components/user/models';
import { UserContextProvider } from '@/components/user/context';


const USER: User = {
  name: faker.person.fullName(),
  email: faker.internet.email(),
  title: faker.person.jobTitle(),
  role: faker.person.jobType(),
  imageUrl: faker.helpers.maybe(() => faker.image.avatar()),
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
