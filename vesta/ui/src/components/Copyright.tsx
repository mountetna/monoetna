import * as React from 'react';
import Typography from '@mui/material/Typography';
import Link from 'next/link';

export default function Copyright() {
  return (
    <Typography variant="pBoldMonoCaps" color="green.grade75" align="center">
      {'Copyright © '}
      <Link color="inherit" href="https://mui.com/">
        Your Website
      </Link>{' '}
      {new Date().getFullYear()}.
    </Typography>
  );
}
