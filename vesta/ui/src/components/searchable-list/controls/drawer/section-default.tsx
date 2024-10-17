'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import { useTheme } from '@mui/material';
import _ from 'lodash'

import { DrawerSectionProps } from './models';
import DrawerPill from './pill';


export default function DrawerSectionDefault({
    name,
    items,
    activeKeys,
    onClickItem,
}: DrawerSectionProps) {
    const theme = useTheme()

    return (
        <Box
            sx={{
                display: 'flex',
                flexDirection: 'column',
                gap: '16px',
                pr: '10px',
            }}
        >
            <Typography
                className='drawer-section-name'
                variant='pMediumBoldWt'
            >
                {_.startCase(name)}
            </Typography>

            <Box
                sx={{
                    display: 'flex',
                    flexDirection: 'row',
                    flexWrap: 'wrap',
                    gap: '10px',
                }}
            >
                {items.map(item => (
                    <DrawerPill
                        variant='yellow'
                        key={item.key}
                        label={item.label}
                        active={activeKeys.has(item.key)}
                        onClick={() => onClickItem(item)}
                    />
                ))}
            </Box>
        </Box>
    )
}