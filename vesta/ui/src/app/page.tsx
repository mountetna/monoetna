import * as React from 'react';
import Container from '@mui/system/Container';
import Typography from '@mui/material/Typography';
import Box from '@mui/system/Box';
import NextLink from 'next/link';
// import ProTip from '@/components/ProTip';
// import Copyright from '@/components/Copyright';

export default function Home() {
  return (
    <Container maxWidth="desktopLg">
      <Box
        sx={{
          my: 4,
          display: 'flex',
          flexDirection: 'column',
          justifyContent: 'center',
          alignItems: 'center',
        }}
      >
        <Typography variant="h1" sx={{ mb: 2 }}>
          Material UI - Next.js App Router example in TypeScript
        </Typography>
        <NextLink href="/about" color="secondary">
          Go to the about page
        </NextLink>
      </Box>
    </Container>
  );
}
