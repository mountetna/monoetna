'use client'

import * as React from 'react';
import { Fade, Tooltip as MUITooltip, styled, tooltipClasses, TooltipProps } from '@mui/material';
import theme from '@/theme';


export default styled(({ className, ...props }: TooltipProps) => (
    <MUITooltip
        {...props}
        TransitionComponent={Fade}
        TransitionProps={{
            easing: theme.transitions.easing.ease,
            timeout: theme.transitions.duration.ease,
        }}
        classes={{ popper: className }}
    />
))(() => ({
    [`& .${tooltipClasses.tooltip}`]: {
        m: '0px',
        p: '0px',
        backgroundColor: 'transparent',
        // maxWidth: 'none,
        border: 'none',
    },
}));