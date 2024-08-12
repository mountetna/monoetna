'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { SxProps, useTheme } from '@mui/material';
import Image, { StaticImageData } from 'next/image';

import { TypographyVariant } from '@/lib/utils/types';
import { useWindowDimensions } from '@/lib/utils/responsive';


export enum Classes {
    base = 'pill',
    filled = 'pill-filled',
    stroked = 'pill-stroked',
}


export default function Pill({
    label,
    typographyVariant = 'pMedium',
    icon,
    iconAlt,
    iconPosition = 'before',
    variant,
    sx = {},
}: {
    label: string,
    typographyVariant: TypographyVariant,
    icon?: StaticImageData,
    iconAlt?: string,
    iconPosition?: 'before' | 'after',
    variant: 'filled' | 'stroked',
    sx?: SxProps,
}) {
    const theme = useTheme()

    // Manage width
    const [isTextWrapped, setIsTextWrapped] = React.useState(false)
    const testLabelRef = React.useRef<HTMLElement>()
    const measureLineHeightRef = React.useRef<HTMLElement>()

    // const {isResizing: isWindowResizing} = useWindowDimensions()

    // React.useEffect(() => {
    //     const testLabelEl = testLabelRef.current
    //     const measureLineHeightEl = measureLineHeightRef.current

    //     if (testLabelEl && measureLineHeightEl &&
    //         testLabelEl.offsetHeight > measureLineHeightEl.offsetHeight
    //     ) {
    //         setIsTextWrapped(true)
    //     } else {
    //         setIsTextWrapped(false)
    //     }
    // }, [isWindowResizing])

    let iconEl
    if (icon) {
        iconEl = (
            <Image
                src={icon}
                alt={iconAlt ? iconAlt : `Image for ${label}`}
            />
        )
    }

    const classes = [Classes.base]
    switch (variant) {
        case 'filled':
            classes.push(Classes.filled)
            break
        case 'stroked':
            classes.push(Classes.stroked)
            break
    }

    const labelEl = (
        <Typography
            variant={typographyVariant}
            sx={{
                display: 'inline-flex',

            }}
        >
            {label}
        </Typography>
    )

    return (
        <Box
            className={classes.join(' ')}
            sx={{
                position: 'relative',
                display: 'flex',
                flexDirection: 'row',
                gap: '6px',
                px: '10px',
                py: '2px',
                borderRadius: '20px',
                alignItems: 'center',
                textAlign: 'center',
                [`&.${Classes.filled}`]: {
                    background: theme.palette.utilityWhite.main,
                },
                [`&.${Classes.stroked}`]: {
                    border: `1px solid ${theme.palette.ground.grade10}`,
                    background: 'transparent',
                },
                ...sx,
            }}
        >
            {iconPosition === 'before' && iconEl}

            <Box
                sx={{
                    '& > *': {
                        width: isTextWrapped ? 'min-content' : 'unset',
                    },
                }}
            >
                {labelEl}
            </Box>

            <Box
                ref={testLabelRef}
                sx={{
                    position: 'absolute',
                    visibility: 'hidden',
                }}
            >
                {labelEl}
            </Box>

            <Box
                ref={measureLineHeightRef}
                sx={{
                    position: 'absolute',
                    visibility: 'hidden',
                }}
            >
                <Typography
                    variant={typographyVariant}
                    sx={{
                        display: 'inline-flex',
                    }}
                >
                    <br />
                </Typography>
            </Box>

            {iconPosition === 'after' && iconEl}
        </Box>
    )
}