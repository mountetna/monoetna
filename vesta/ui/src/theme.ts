'use client';
import { createTheme } from '@mui/material/styles';
import { theFuture, relativeMonoPro10Pitch, relativeMonoPro11Pitch, relativeMonoPro12Pitch } from '@/fonts';


// https://mui.com/material-ui/customization/typography/#adding-amp-disabling-variants
declare module '@mui/material/styles' {
  interface TypographyVariants {
    pTitleRegularWt: React.CSSProperties;
    pTitleMediumWt: React.CSSProperties;
    pTitleMono: React.CSSProperties;
    pTitleMonoCaps: React.CSSProperties;
    pSubtitle: React.CSSProperties;
    pSubtitleBoldWt: React.CSSProperties;
    pSubtitleMono: React.CSSProperties;
    pSubtitleMonoCaps: React.CSSProperties;
    pMedium: React.CSSProperties;
    pMediumMediumWt: React.CSSProperties;
    pMediumBoltWt: React.CSSProperties;
    pMediumMono: React.CSSProperties;
    pMediumMonoCaps: React.CSSProperties;
    pBoldMonoCaps: React.CSSProperties;
    pBody: React.CSSProperties;
    pBodyMediumWt: React.CSSProperties;
    pBodyBoldWt: React.CSSProperties;
    pBodyMono: React.CSSProperties;
    pXS: React.CSSProperties;
    p2XS: React.CSSProperties;
    p2XSBoldWt: React.CSSProperties;
    p3XS: React.CSSProperties;
    p3XSMediumWt: React.CSSProperties;
    p3XSBoldWt: React.CSSProperties;
    p3XSMono: React.CSSProperties;
    p4XS: React.CSSProperties;
    pLarge: React.CSSProperties;
    pLargeMediumWt: React.CSSProperties;
    pLargeBoldWt: React.CSSProperties;
  }

  // allow configuration using `createTheme`
  interface TypographyVariantsOptions {
    pTitleRegularWt?: React.CSSProperties;
    pTitleMediumWt?: React.CSSProperties;
    pTitleMono?: React.CSSProperties;
    pTitleMonoCaps?: React.CSSProperties;
    pSubtitle?: React.CSSProperties;
    pSubtitleBoldWt?: React.CSSProperties;
    pSubtitleMono?: React.CSSProperties;
    pSubtitleMonoCaps?: React.CSSProperties;
    pMedium?: React.CSSProperties;
    pMediumMediumWt?: React.CSSProperties;
    pMediumBoltWt?: React.CSSProperties;
    pMediumMono?: React.CSSProperties;
    pMediumMonoCaps?: React.CSSProperties;
    pBoldMonoCaps?: React.CSSProperties;
    pBody?: React.CSSProperties;
    pBodyMediumWt?: React.CSSProperties;
    pBodyBoldWt?: React.CSSProperties;
    pBodyMono?: React.CSSProperties;
    pXS?: React.CSSProperties;
    p2XS?: React.CSSProperties;
    p2XSBoldWt?: React.CSSProperties;
    p3XS?: React.CSSProperties;
    p3XSMediumWt?: React.CSSProperties;
    p3XSBoldWt?: React.CSSProperties;
    p3XSMono?: React.CSSProperties;
    p4XS?: React.CSSProperties;
    pLarge?: React.CSSProperties;
    pLargeMediumWt?: React.CSSProperties;
    pLargeBoldWt?: React.CSSProperties;
  }

  interface BreakpointOverrides {
    mobile: true;
    tablet: true;
    desktop: true;
    desktopLg: true;
    xs: false;
    sm: false;
    md: false;
    lg: false;
    xl: false;
  }
}

// Update the Typography's variant prop options
declare module '@mui/material/Typography' {
  interface TypographyPropsVariantOverrides {
    pTitleRegularWt: true;
    pTitleMediumWt: true;
    pTitleMono: true;
    pTitleMonoCaps: true;
    pSubtitle: true;
    pSubtitleBoldWt: true;
    pSubtitleMono: true;
    pSubtitleMonoCaps: true;
    pMedium: true;
    pMediumMediumWt: true;
    pMediumBoltWt: true;
    pMediumMono: true;
    pMediumMonoCaps: true;
    pBoldMonoCaps: true;
    pBody: true;
    pBodyMediumWt: true;
    pBodyBoldWt: true;
    pBodyMono: true;
    pXS: true;
    p2XS: true;
    p2XSBoldWt: true;
    p3XS: true;
    p3XSMediumWt: true;
    p3XSBoldWt: true;
    p3XSMono: true;
    p4XS: true;
    pLarge: true;
    pLargeMediumWt: true;
    pLargeBoldWt: true;
    subtitle1: false;
    subtitle2: false;
    body1: false;
    body2: false;
    button: false;
    caption: false;
    overline: false;
  }
}


const LETTER_SPACING = '-0.01em'

const theme = createTheme({
  palette: {
    mode: 'light',
  },
  breakpoints: {
    values: {
      mobile: 0,
      tablet: 700,
      desktop: 1200,
      desktopLg: 1440,
    }
  },
  typography: {
    // https://mui.com/material-ui/customization/typography/#variants
    fontFamily: theFuture.style.fontFamily,
    // h1-6, pLarge defined below
    pTitleRegularWt: {
      fontSize: 22,
      lineHeight: '140%',
      letterSpacing: LETTER_SPACING,
    },
    pTitleMediumWt: {
      fontSize: 22,
      fontWeight: 500,
      lineHeight: '140%',
      letterSpacing: LETTER_SPACING,
    },
    pTitleMono: {
      fontFamily: relativeMonoPro10Pitch.style.fontFamily,
      fontSize: 22,
      lineHeight: '140%',
      letterSpacing: LETTER_SPACING,
    },
    pTitleMonoCaps: {
      fontFamily: relativeMonoPro10Pitch.style.fontFamily,
      fontSize: 22,
      lineHeight: '140%',
      letterSpacing: LETTER_SPACING,
      textTransform: 'uppercase',
    },
    pSubtitle: {
      fontSize: 20,
      fontWeight: 500,
      lineHeight: '140%',
      letterSpacing: LETTER_SPACING,
    },
    pSubtitleBoldWt: {
      fontSize: 20,
      fontWeight: 700,
      lineHeight: '140%',
      letterSpacing: LETTER_SPACING,
    },
    pSubtitleMono: {
      fontFamily: relativeMonoPro10Pitch.style.fontFamily,
      fontSize: 20,
      lineHeight: '140%',
      letterSpacing: LETTER_SPACING,
    },
    pSubtitleMonoCaps: {
      fontFamily: relativeMonoPro10Pitch.style.fontFamily,
      fontSize: 20,
      lineHeight: '140%',
      letterSpacing: LETTER_SPACING,
      textTransform: 'uppercase',
    },
    pMedium: {
      fontSize: 18,
      lineHeight: '150%',
      letterSpacing: LETTER_SPACING,
    },
    pMediumMediumWt: {
      fontSize: 18,
      fontWeight: 500,
      lineHeight: '150%',
      letterSpacing: LETTER_SPACING,
    },
    pMediumBoltWt: {
      fontSize: 18,
      fontWeight: 700,
      lineHeight: '150%',
      letterSpacing: LETTER_SPACING,
    },
    pMediumMono: {
      fontFamily: relativeMonoPro10Pitch.style.fontFamily,
      fontSize: 18,
      lineHeight: '150%',
      letterSpacing: LETTER_SPACING,
    },
    pMediumMonoCaps: {
      fontFamily: relativeMonoPro10Pitch.style.fontFamily,
      fontSize: 18,
      fontWeight: 500,
      lineHeight: '150%',
      letterSpacing: LETTER_SPACING,
      textTransform: 'uppercase',
    },
    pBoldMonoCaps: {
      fontFamily: relativeMonoPro10Pitch.style.fontFamily,
      fontSize: 18,
      fontWeight: 700,
      lineHeight: '153%',
      letterSpacing: LETTER_SPACING,
      textTransform: 'uppercase',
    },
    pBody: {
      fontSize: 16,
      lineHeight: '153%',
      letterSpacing: LETTER_SPACING,
    },
    pBodyMediumWt: {
      fontSize: 16,
      fontWeight: 500,
      lineHeight: '153%',
      letterSpacing: LETTER_SPACING,
    },
    pBodyBoldWt: {
      fontSize: 16,
      fontWeight: 700,
      lineHeight: '153%',
      letterSpacing: LETTER_SPACING,
    },
    pBodyMono: {
      fontFamily: relativeMonoPro10Pitch.style.fontFamily,
      fontSize: 16,
      lineHeight: '153%',
      letterSpacing: LETTER_SPACING,
    },
    pXS: {
      fontSize: 15,
      lineHeight: '153%',
      letterSpacing: LETTER_SPACING,
    },
    p2XS: {
      fontSize: 14,
      lineHeight: '153%',
      letterSpacing: LETTER_SPACING,
    },
    p2XSBoldWt: {
      fontSize: 14,
      fontWeight: 700,
      lineHeight: '153%',
      letterSpacing: LETTER_SPACING,
    },
    p3XS: {
      fontSize: 12,
      lineHeight: '122%',
      letterSpacing: LETTER_SPACING,
    },
    p3XSMediumWt: {
      fontSize: 12,
      fontWeight: 500,
      lineHeight: '122%',
      letterSpacing: LETTER_SPACING,
    },
    p3XSBoldWt: {
      fontSize: 12,
      fontWeight: 700,
      lineHeight: '122%',
      letterSpacing: LETTER_SPACING,
    },
    p3XSMono: {
      fontFamily: relativeMonoPro10Pitch.style.fontFamily,
      fontSize: 12,
      fontWeight: 700,
      lineHeight: '122%',
      letterSpacing: LETTER_SPACING,
    },
    p4XS: {
      fontSize: 10,
      lineHeight: '122%',
      letterSpacing: LETTER_SPACING,
    },
    subtitle1: undefined,
    subtitle2: undefined,
    body1: undefined,
    body2: undefined,
    button: undefined,
    caption: undefined,
    overline: undefined,
  },
  components: {},
});

theme.typography.pLarge = {
  fontFamily: theFuture.style.fontFamily,
  fontSize: 18,
  lineHeight: '150%',
  letterSpacing: LETTER_SPACING,
  [theme.breakpoints.up('tablet')]: {
    fontSize: 20,
    lineHeight: '148%',
  },
  [theme.breakpoints.up('desktop')]: {
    fontSize: 22,
    lineHeight: '148%',
  },
}

// only for tablet, desktop
theme.typography.pLargeMediumWt = {
  fontFamily: theFuture.style.fontFamily,
  fontSize: 20,
  fontWeight: 500,
  lineHeight: '148%',
  letterSpacing: LETTER_SPACING,
  [theme.breakpoints.up('desktop')]: {
    fontSize: 22,
    lineHeight: '148%',
  },
}

// only for mobile, tablet
theme.typography.pLargeBoldWt = {
  fontFamily: theFuture.style.fontFamily,
  fontSize: 18,
  fontWeight: 700,
  lineHeight: '150%',
  letterSpacing: LETTER_SPACING,
  [theme.breakpoints.up('desktop')]: {
    fontSize: 20,
    lineHeight: '148%',
  },
}

theme.typography.h1 = {
  fontFamily: theFuture.style.fontFamily,
  fontSize: 36,
  fontWeight: 700,
  lineHeight: '135%',
  letterSpacing: LETTER_SPACING,
  [theme.breakpoints.up('tablet')]: {
    fontSize: 46,
    lineHeight: '140%',
  },
  [theme.breakpoints.up('desktop')]: {
    fontSize: 56,
    lineHeight: '144%',
  },
}

theme.typography.h2 = {
  fontFamily: theFuture.style.fontFamily,
  fontSize: 24,
  fontWeight: 700,
  lineHeight: '140%',
  letterSpacing: LETTER_SPACING,
  [theme.breakpoints.up('tablet')]: {
    fontSize: 36,
    lineHeight: '135%',
  },
  [theme.breakpoints.up('desktop')]: {
    fontSize: 48,
    lineHeight: '144%',
  },
}

export default theme;
