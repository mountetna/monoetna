'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import ButtonBase from '@mui/material/ButtonBase';
import Collapse from '@mui/material/Collapse';
import Fade from '@mui/material/Fade';
import { useTheme } from '@mui/material';
import { TransitionProps } from '@mui/material/transitions';
import _ from 'lodash'
import Image from 'next/image';

import { DrawerSectionProps } from './models';

import indicatorArrowDark from '/public/images/icons/indicator-arrow-dark.svg'
import DrawerPill from './pill';


const activeClass = 'active'

interface Props extends DrawerSectionProps {
    open: boolean;
    onSetOpen: (open: boolean) => void;
}


export default function DrawerSectionExpandable({
    name,
    items,
    activeKeys,
    onClickItem,
    open,
    onSetOpen,
}: Props) {
    const theme = useTheme()

    const animationProps: TransitionProps = {
        in: open,
        easing: theme.transitions.easing.quint,
        timeout: theme.transitions.duration.quint,
    }

    return (
        <Box
            sx={{
                p: '2px',
            }}
        >
            <ButtonBase
                onClick={() => onSetOpen(!open)}
                sx={{
                    display: 'flex',
                    justifyContent: 'space-between',
                    width: '100%',
                }}
            >
                <Typography
                    className='drawer-section-name'
                    variant='pMediumBoldWt'
                    sx={{
                        color: 'ground.grade10',
                    }}
                >
                    {_.startCase(name)}
                </Typography>

                <Image
                    src={indicatorArrowDark}
                    alt={`Triangle icon pointing ${open ? 'up' : 'down'}`}
                    style={{
                        rotate: `${open ? 0 : -180}deg`,
                        transition: theme.transitions.create(
                            'all',
                            {
                                easing: theme.transitions.easing.quint,
                                duration: theme.transitions.duration.quint,
                            },
                        ),
                    }}
                />
            </ButtonBase>

            <Collapse
                {...animationProps}
            >
                <Fade
                    {...animationProps}
                >
                    <Box
                        sx={{
                            display: 'flex',
                            flexDirection: 'row',
                            flexWrap: 'wrap',
                            gap: '10px',
                            pt: '16px',
                        }}
                    >
                        {items.map(item => (
                            <DrawerPill
                                key={item.key}
                                label={item.label}
                                active={activeKeys.has(item.key)}
                                onClick={() => onClickItem(item)}
                            />
                        ))}
                    </Box>
                </Fade>
            </Collapse>
        </Box>
    )
}