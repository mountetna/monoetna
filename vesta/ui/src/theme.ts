'use client';

import { createTheme, alpha, PaletteColor, SxProps, Theme } from '@mui/material/styles';
// import { PaletteColorOptions } from '@mui/material/styles/createPalette';
import { theFuture, relativeMonoPro10Pitch } from '@/fonts';


declare module '@mui/material/styles' {
  // https://mui.com/material-ui/customization/palette/#typescript
  interface PalleteColorsDL {
    grade100: string;
    grade75: string;
    grade50: string;
    grade25: string;
    grade10: string;
  }

  interface Palette {
    utilityLowlight: Palette['primary'];
    utilityHighlight: Palette['primary'];
    utilityWhite: Palette['primary'];
    utilityWhiteTransparent25: Palette['primary'];
    utilityUCSFNavy: Palette['primary'];
    utilityUCSFLightNavy: Palette['primary'];
    ground: Palette['primary'] & PalleteColorsDL;
    green: Palette['primary'] & PalleteColorsDL;
    teal: Palette['primary'] & PalleteColorsDL;
    blue: Palette['primary'] & PalleteColorsDL;
    magenta: Palette['primary'] & PalleteColorsDL;
    red: Palette['primary'] & PalleteColorsDL;
    orange: Palette['primary'] & PalleteColorsDL;
    yellow: Palette['primary'] & PalleteColorsDL;
  }

  interface PalleteColorOptionsDL {
    grade100?: string;
    grade75?: string;
    grade50?: string;
    grade25?: string;
    grade10?: string;
  }

  interface PaletteOptions {
    utilityLowlight?: PaletteOptions['primary'];
    utilityHighlight?: PaletteOptions['primary'];
    utilityWhite?: PaletteOptions['primary'];
    utilityWhiteTransparent25?: PaletteOptions['primary'];
    utilityUCSFNavy?: PaletteOptions['primary'];
    utilityUCSFLightNavy?: PaletteOptions['primary'];
    ground?: PaletteOptions['primary'] & PalleteColorOptionsDL;
    green?: PaletteOptions['primary'] & PalleteColorOptionsDL;
    teal?: PaletteOptions['primary'] & PalleteColorOptionsDL;
    blue?: PaletteOptions['primary'] & PalleteColorOptionsDL;
    magenta?: PaletteOptions['primary'] & PalleteColorOptionsDL;
    red?: PaletteOptions['primary'] & PalleteColorOptionsDL;
    orange?: PaletteOptions['primary'] & PalleteColorOptionsDL;
    yellow?: PaletteOptions['primary'] & PalleteColorOptionsDL;
  }

  // https://mui.com/material-ui/customization/typography/#adding-amp-disabling-variants
  interface TypographyVariants {
    h3Digits: React.CSSProperties;
    h3MediumWt: React.CSSProperties;
    h3SmallCaps: React.CSSProperties;
    h4SmallCaps: React.CSSProperties;
    h5BoldWt: React.CSSProperties;
    h6BoldWt: React.CSSProperties;
    h6SmallCaps: React.CSSProperties;
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
    pMediumBoldWt: React.CSSProperties;
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
    h3Digits?: React.CSSProperties;
    h3MediumWt?: React.CSSProperties;
    h3SmallCaps?: React.CSSProperties;
    h4SmallCaps?: React.CSSProperties;
    h5BoldWt?: React.CSSProperties;
    h6BoldWt?: React.CSSProperties;
    h6SmallCaps?: React.CSSProperties;
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
    pMediumBoldWt?: React.CSSProperties;
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

  interface Duration {
    ease: number;
    swell: number;
    expo: number;
    quint: number;
    long: number;
  }

  interface Easing {
    ease: string;
    swell: string;
    expo: string;
    quint: string;
  }
}

// Update the Typography's variant prop options
declare module '@mui/material/Typography' {
  interface TypographyPropsVariantOverrides {
    h3Digits: true;
    h3MediumWt: true;
    h3SmallCaps: true;
    h4SmallCaps: true;
    h5BoldWt: true;
    h6BoldWt: true;
    h6SmallCaps: true;
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
    pMediumBoldWt: true;
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

const DEFAULT_LETTER_SPACING = '-0.01em'
const DEFAULT_CONTRAST_RATIO = 4.5

let theme = createTheme({
  palette: {
    mode: 'light',
    // https://mui.com/material-ui/customization/palette/#accessibility
    contrastThreshold: DEFAULT_CONTRAST_RATIO,
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
      fontFamily: theFuture.style.fontFamily,
      fontSize: 22,
      lineHeight: '140%',
      letterSpacing: DEFAULT_LETTER_SPACING,
    },
    pTitleMediumWt: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 22,
      fontWeight: 500,
      lineHeight: '140%',
      letterSpacing: DEFAULT_LETTER_SPACING,
    },
    pTitleMono: {
      fontFamily: relativeMonoPro10Pitch.style.fontFamily,
      fontSize: 22,
      lineHeight: '140%',
      letterSpacing: DEFAULT_LETTER_SPACING,
    },
    pTitleMonoCaps: {
      fontFamily: relativeMonoPro10Pitch.style.fontFamily,
      fontSize: 22,
      lineHeight: '140%',
      letterSpacing: DEFAULT_LETTER_SPACING,
      textTransform: 'uppercase',
    },
    pSubtitle: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 20,
      fontWeight: 500,
      lineHeight: '140%',
      letterSpacing: DEFAULT_LETTER_SPACING,
    },
    pSubtitleBoldWt: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 20,
      fontWeight: 700,
      lineHeight: '140%',
      letterSpacing: DEFAULT_LETTER_SPACING,
    },
    pSubtitleMono: {
      fontFamily: relativeMonoPro10Pitch.style.fontFamily,
      fontSize: 20,
      lineHeight: '140%',
      letterSpacing: DEFAULT_LETTER_SPACING,
    },
    pSubtitleMonoCaps: {
      fontFamily: relativeMonoPro10Pitch.style.fontFamily,
      fontSize: 20,
      lineHeight: '140%',
      letterSpacing: DEFAULT_LETTER_SPACING,
      textTransform: 'uppercase',
    },
    pMedium: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 18,
      lineHeight: '150%',
      letterSpacing: DEFAULT_LETTER_SPACING,
    },
    pMediumMediumWt: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 18,
      fontWeight: 500,
      lineHeight: '150%',
      letterSpacing: DEFAULT_LETTER_SPACING,
    },
    pMediumBoldWt: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 18,
      fontWeight: 700,
      lineHeight: '150%',
      letterSpacing: DEFAULT_LETTER_SPACING,
    },
    pMediumMono: {
      fontFamily: relativeMonoPro10Pitch.style.fontFamily,
      fontSize: 18,
      lineHeight: '150%',
      letterSpacing: DEFAULT_LETTER_SPACING,
    },
    pMediumMonoCaps: {
      fontFamily: relativeMonoPro10Pitch.style.fontFamily,
      fontSize: 18,
      fontWeight: 500,
      lineHeight: '150%',
      letterSpacing: DEFAULT_LETTER_SPACING,
      textTransform: 'uppercase',
    },
    pBoldMonoCaps: {
      fontFamily: relativeMonoPro10Pitch.style.fontFamily,
      fontSize: 18,
      fontWeight: 700,
      lineHeight: '153%',
      letterSpacing: DEFAULT_LETTER_SPACING,
      textTransform: 'uppercase',
    },
    pBody: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 16,
      lineHeight: '153%',
      letterSpacing: DEFAULT_LETTER_SPACING,
    },
    pBodyMediumWt: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 16,
      fontWeight: 500,
      lineHeight: '153%',
      letterSpacing: DEFAULT_LETTER_SPACING,
    },
    pBodyBoldWt: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 16,
      fontWeight: 700,
      lineHeight: '153%',
      letterSpacing: DEFAULT_LETTER_SPACING,
    },
    pBodyMono: {
      fontFamily: relativeMonoPro10Pitch.style.fontFamily,
      fontSize: 16,
      lineHeight: '153%',
      letterSpacing: DEFAULT_LETTER_SPACING,
    },
    pXS: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 15,
      lineHeight: '153%',
      letterSpacing: DEFAULT_LETTER_SPACING,
    },
    p2XS: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 14,
      lineHeight: '153%',
      letterSpacing: DEFAULT_LETTER_SPACING,
    },
    p2XSBoldWt: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 14,
      fontWeight: 700,
      lineHeight: '153%',
      letterSpacing: DEFAULT_LETTER_SPACING,
    },
    p3XS: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 12,
      lineHeight: '122%',
      letterSpacing: DEFAULT_LETTER_SPACING,
    },
    p3XSMediumWt: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 12,
      fontWeight: 500,
      lineHeight: '122%',
      letterSpacing: DEFAULT_LETTER_SPACING,
    },
    p3XSBoldWt: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 12,
      fontWeight: 700,
      lineHeight: '122%',
      letterSpacing: DEFAULT_LETTER_SPACING,
    },
    p3XSMono: {
      fontFamily: relativeMonoPro10Pitch.style.fontFamily,
      fontSize: 12,
      fontWeight: 700,
      lineHeight: '122%',
      letterSpacing: DEFAULT_LETTER_SPACING,
    },
    p4XS: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 10,
      lineHeight: '122%',
      letterSpacing: DEFAULT_LETTER_SPACING,
    },
    subtitle1: undefined,
    subtitle2: undefined,
    body1: undefined,
    body2: undefined,
    button: undefined,
    caption: undefined,
    overline: undefined,
  },
  transitions: {
    duration: {
      ease: 300,
      swell: 350,
      expo: 360,
      quint: 450,
      long: 2000,
    },
    easing: {
      ease: 'ease',
      swell: 'cubic-bezier(0.85, 0, 0.15, 1)',
      expo: 'cubic-bezier(0.16, 1, 0.3, 1)',
      quint: 'cubic-bezier(0.92, 0, 0.06, 1)',
    },
  },
  components: {
    // https://mui.com/material-ui/customization/theme-components/#theme-default-props
    MuiContainer: {
      defaultProps: {
        maxWidth: 'desktopLg',
      },
    },
    MuiButtonBase: {
      defaultProps: {
        disableRipple: true,
      },
    },
  }
});

interface ColorTokens {
  main: string
  [key: string]: string
}

const colors: Record<string, string | ColorTokens> = {
  utilityLowlight: "#121212",
  utilityHighlight: '#F9F8F6',
  utilityWhite: '#FFFFFF',
  utilityWhiteTransparent25: alpha('#FFFFFF', 0.75),
  utilityUCSFNavy: '#052049',
  utilityUCSFLightNavy: '#E6F4FB',
  ground: {
    main: '#9F9F9C',
    grade100: '#F3F2E3',
    grade75: '#DDD9D0',
    grade50: '#9F9F9C',
    grade25: '#3B3B3B',
    grade10: '#242424',
  },
  green: {
    main: '#6C6D2A',
    grade100: '#D6D8A8',
    grade75: '#A2A648',
    grade50: '#6C6D2A',
    grade25: '#3D4A1C',
    grade10: '#253413',
  },
  teal: {
    main: '#556E66',
    grade100: '#B8CBC2',
    grade75: '#7FA190',
    grade50: '#556E66',
    grade25: '#3A4B48',
    grade10: '#283234',
  },
  blue: {
    main: '#519BCB',
    grade100: '#C2D0E4',
    grade75: '#89A7CE',
    grade50: '#519BCB',
    grade25: '#2A5A8D',
    grade10: '#13283F',
  },
  magenta: {
    main: '#E793B0',
    grade100: '#EEDFE4',
    grade75: '#E4B8C7',
    grade50: '#E793B0',
    grade25: '#C46485',
    grade10: '#4A1C1C',
  },
  red: {
    main: '#DB7976',
    grade100: '#EAD7D3',
    grade75: '#D6AFA6',
    grade50: '#DB7976',
    grade25: '#B53B38',
    grade10: '#4B1716',
  },
  orange: {
    main: '#D36F49',
    grade100: '#F0D7C2',
    grade75: '#DDA373',
    grade50: '#D36F49',
    grade25: '#9C4626',
    grade10: '#562915',
  },
  yellow: {
    main: '#E9C54E',
    grade100: '#FBEEC6',
    grade75: '#F7DE8E',
    grade50: '#E9C54E',
    grade25: '#AA832C',
    grade10: '#3E3528',
  },
}

// https://mui.com/material-ui/customization/palette/#generate-tokens-using-augmentcolor-utility
const paletteColors: Record<string, PaletteColor> = {}
for (const [name, val] of Object.entries(colors)) {

  let colorTokens: ColorTokens;
  if (typeof val == 'string') {
    colorTokens = { main: val }
  } else {
    colorTokens = val
  }

  // generate standard tokens
  paletteColors[name] = theme.palette.augmentColor({
    color: {
      main: colorTokens.main,
    },
    name: name,
  })

  // add custom tokens
  paletteColors[name] = { ...paletteColors[name], ...colorTokens }
}

export const containerPadding = {
  padding: '0 8px',
  [theme.breakpoints.up('tablet')]: {
    padding: '0 16px',
  },
  [theme.breakpoints.up('desktop')]: {
    padding: '0 16px',
  },
}

theme = createTheme(theme, {
  palette: {
    ...paletteColors,
    background: {
      default: paletteColors.utilityHighlight.main,
    },
  },
  typography: {
    h1: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 36,
      fontWeight: 700,
      lineHeight: '135%',
      letterSpacing: DEFAULT_LETTER_SPACING,
      [theme.breakpoints.up('tablet')]: {
        fontSize: 46,
        lineHeight: '140%',
      },
      [theme.breakpoints.up('desktop')]: {
        fontSize: 56,
        lineHeight: '144%',
      },
    },
    h2: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 24,
      fontWeight: 700,
      lineHeight: '140%',
      letterSpacing: DEFAULT_LETTER_SPACING,
      [theme.breakpoints.up('tablet')]: {
        fontSize: 36,
        lineHeight: '135%',
      },
      [theme.breakpoints.up('desktop')]: {
        fontSize: 48,
        lineHeight: '144%',
      },
    },
    h3: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 22,
      fontWeight: 700,
      lineHeight: '140%',
      letterSpacing: DEFAULT_LETTER_SPACING,
      [theme.breakpoints.up('tablet')]: {
        fontSize: 36,
        lineHeight: '140%',
      },
      [theme.breakpoints.up('desktop')]: {
        fontSize: 46,
        lineHeight: '140%',
      },
    },
    h3MediumWt: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 22,
      fontWeight: 500,
      lineHeight: '140%',
      letterSpacing: DEFAULT_LETTER_SPACING,
      [theme.breakpoints.up('tablet')]: {
        fontSize: 36,
        lineHeight: '140%',
      },
      [theme.breakpoints.up('desktop')]: {
        fontSize: 46,
        lineHeight: '140%',
      },
    },
    h3SmallCaps: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 22,
      fontWeight: 700,
      lineHeight: '140%',
      letterSpacing: DEFAULT_LETTER_SPACING,
      fontVariant: 'all-small-caps', // TODO: verify this is accurate
    },
    h3Digits: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 46,
      fontWeight: 700,
      lineHeight: '140%',
      letterSpacing: DEFAULT_LETTER_SPACING,
      [theme.breakpoints.up('tablet')]: {
        fontSize: 46,
        lineHeight: '140%',
      },
      [theme.breakpoints.up('desktop')]: {
        fontSize: 46,
        lineHeight: '140%',
      },
    },
    h4: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 22,
      fontWeight: 700,
      lineHeight: '140%',
      letterSpacing: DEFAULT_LETTER_SPACING,
      [theme.breakpoints.up('tablet')]: {
        fontSize: 24,
        lineHeight: '140%',
      },
      [theme.breakpoints.up('desktop')]: {
        fontSize: 36,
        lineHeight: '135%',
      },
    },
    h4SmallCaps: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 18,
      fontWeight: 700,
      lineHeight: '135%',
      letterSpacing: '0.04em',
      fontVariant: 'all-small-caps', // TODO: verify this is accurate
    },
    h5: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 20,
      fontWeight: 500,
      lineHeight: '140%',
      letterSpacing: DEFAULT_LETTER_SPACING,
      [theme.breakpoints.up('tablet')]: {
        fontSize: 22,
        lineHeight: '140%',
      },
      [theme.breakpoints.up('desktop')]: {
        fontSize: 24,
        lineHeight: '135%',
      },
    },
    h5BoldWt: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 20,
      fontWeight: 700,
      lineHeight: '140%',
      letterSpacing: DEFAULT_LETTER_SPACING,
      [theme.breakpoints.up('tablet')]: {
        fontSize: 22,
        lineHeight: '140%',
      },
      [theme.breakpoints.up('desktop')]: {
        fontSize: 24,
        lineHeight: '135%',
      },
    },
    h6: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 16,
      fontWeight: 500,
      lineHeight: '140%',
      letterSpacing: DEFAULT_LETTER_SPACING,
      [theme.breakpoints.up('tablet')]: {
        fontSize: 18,
        lineHeight: '140%',
      },
      [theme.breakpoints.up('desktop')]: {
        fontSize: 22,
        lineHeight: '135%',
      },
    },
    h6BoldWt: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 18,
      fontWeight: 700,
      lineHeight: '140%',
      letterSpacing: DEFAULT_LETTER_SPACING,
      [theme.breakpoints.up('tablet')]: {
        fontSize: 18,
        lineHeight: '140%',
      },
      [theme.breakpoints.up('desktop')]: {
        fontSize: 22,
        lineHeight: '135%',
      },
    },
    h6SmallCaps: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 14,
      fontWeight: 700,
      lineHeight: '153%',
      letterSpacing: '0.05em',
      fontVariant: 'all-small-caps', // TODO: verify this is accurate
    },
    pLarge: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 18,
      lineHeight: '150%',
      letterSpacing: DEFAULT_LETTER_SPACING,
      [theme.breakpoints.up('tablet')]: {
        fontSize: 20,
        lineHeight: '148%',
      },
      [theme.breakpoints.up('desktop')]: {
        fontSize: 22,
        lineHeight: '148%',
      },
    },
    // only for tablet, desktop
    pLargeMediumWt: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 20,
      fontWeight: 500,
      lineHeight: '148%',
      letterSpacing: DEFAULT_LETTER_SPACING,
      [theme.breakpoints.up('desktop')]: {
        fontSize: 22,
        lineHeight: '148%',
      },
    },
    // only for mobile, tablet
    pLargeBoldWt: {
      fontFamily: theFuture.style.fontFamily,
      fontSize: 18,
      fontWeight: 700,
      lineHeight: '150%',
      letterSpacing: DEFAULT_LETTER_SPACING,
      [theme.breakpoints.up('desktop')]: {
        fontSize: 20,
        lineHeight: '148%',
      },
    },
  },
  components: {
    MuiContainer: {
      styleOverrides: {
        root: {...containerPadding},
      },
    },
  }
})

// figure out a better way to incorporate in theme
export const spacing: SxProps = {
  mt: '45px',
  [theme.breakpoints.up('tablet')]: {
    mt: '80px',
  },
  [theme.breakpoints.up('desktop')]: {
    mt: '100px',
  },
  '& > *': {
    mb: '145px',
    [theme.breakpoints.up('tablet')]: {
      mb: '158px',
    },
    [theme.breakpoints.up('desktop')]: {
      mb: '158px',
    },
  },
}

export default theme;