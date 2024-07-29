import * as React from 'react'
import MUILink from '@mui/material/Link';
import Collapse from '@mui/material/Collapse'
import { TransitionProps } from '@mui/material/transitions';
import Fade from '@mui/material/Fade'
import { useTheme } from '@mui/material';
import Box from '@mui/material/Box';
import Link from 'next/link'
import { useRouter } from 'next/navigation';

import LibraryCardButton from '../library-card/library-card-button';


export default function MobileNav({
    open,
}: {
    open: boolean,
}) {
    const theme = useTheme()

    const animationProps: TransitionProps = {
        in: open,
        easing: theme.transitions.easing.quint,
        timeout: theme.transitions.duration.quint,
    }

    return (
        <Collapse
            component='nav'
            aria-label='Main'
            {...animationProps}
        >
            <Fade
                {...animationProps}
            >
                <Box>placeholder</Box>
            </Fade>
        </Collapse>
    )
}