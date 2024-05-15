'use client';
import { createTheme } from '@mui/material/styles';
import { theFuture, relativeMonoPro10Pitch, relativeMonoPro11Pitch, relativeMonoPro12Pitch } from '@/fonts';

const theme = createTheme({
  palette: {
    mode: 'light',
  },
  typography: {
    fontFamily: theFuture.style.fontFamily,
  },
  components: {
    MuiAlert: {
      styleOverrides: {
        root: ({ ownerState }) => ({
          ...(ownerState.severity === 'info' && {
            backgroundColor: '#60a5fa',
          }),
        }),
      },
    },
  },
});

export default theme;
