import * as React from 'react';
import Typography from '@mui/material/Typography';
import Box from '@mui/system/Box';
import NextLink from 'next/link';
import Container from '@mui/system/Container';

import Hero from '@/components/hero/Hero';


export default function Home() {
  return (
    <React.Fragment>
      <Hero />
    </React.Fragment>
  );
}
