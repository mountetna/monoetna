'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import ButtonBase from '@mui/material/ButtonBase';
import { useTheme } from '@mui/material';
import _ from 'lodash'

import { DrawerSectionProps } from './models';


const activeClass = 'active'

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

            }}
        >
            <Typography
                className='drawer-section'
                variant='pMediumBoldWt'
            >
                {_.startCase(name)}
            </Typography>

            <Box>
                {items.map(item => (
                    <ButtonBase
                        key={item.key}
                        className={`drawer-item${activeKeys.has(item.key) ? ` ${activeClass}` : ''}`}
                        onClick={() => onClickItem(item)}
                        sx={{

                        }}
                    >
                        <Typography variant='pMedium'>
                            {item.label}
                        </Typography>
                    </ButtonBase>
                ))}
            </Box>
        </Box>
    )
}