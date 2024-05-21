import * as React from 'react';
import { AppRouterCacheProvider } from '@mui/material-nextjs/v14-appRouter';
import { ThemeProvider } from '@mui/material/styles';
import CssBaseline from '@mui/material/CssBaseline';

import theme from '@/theme';
import Nav from '@/components/nav/Nav';
import Footer from '@/components/Footer';

export default function RootLayout(props: { children: React.ReactNode }) {
  return (
    <html lang="en">
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
