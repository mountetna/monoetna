import * as React from 'react';
import { AppRouterCacheProvider } from '@mui/material-nextjs/v14-appRouter';
import { ThemeProvider } from '@mui/material/styles';
import CssBaseline from '@mui/material/CssBaseline';
import 'swiper/css';
import 'swiper/css/bundle';

import theme from '@/theme';
import Nav from '@/components/nav/nav';
import Footer from '@/components/footer';

export default function RootLayout(props: { children: React.ReactNode }) {
  return (
    <html lang="en">
      <head>
      </head>
      <body>
        <AppRouterCacheProvider>
          <ThemeProvider theme={theme}>
            <CssBaseline enableColorScheme={true} />

            <Nav />
            <main>{props.children}</main>
            <Footer />

          </ThemeProvider>
        </AppRouterCacheProvider>
      </body>
    </html>
  );
}
