import * as React from 'react';
import { Metadata } from 'next';
import Typography from '@mui/material/Typography';
import Box from '@mui/system/Box';
import NextLink from 'next/link';
import Container from '@mui/system/Container';

import Hero from '@/components/hero/Hero';


export const metadata: Metadata = {
  title: 'UCSF Data Library',
}

export default function Home() {
  return (
    <React.Fragment>
      <Hero />
    </React.Fragment>
  );
}
