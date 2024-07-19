'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import ButtonBase from '@mui/material/ButtonBase';
import { useTheme } from '@mui/material';
import { useSpring, animated, useIsomorphicLayoutEffect } from '@react-spring/web';
import _ from 'lodash'
import Image from 'next/image';

import { useWindowDimensions } from '@/lib/utils/responsive';
import { DrawerSectionProps } from './models';

import indicatorArrowDark from '/public/images/icons/indicator-arrow-dark.svg'


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

    const { isResizing: isWindowResizing } = useWindowDimensions()

    // Manage open/close
    const itemsContainerRef = React.useRef<HTMLElement>()
    const [itemsContainerStyle, itemsContainerApi] = useSpring(() => ({
        height: '0px',
        opacity: 0,
    }), [])

    useIsomorphicLayoutEffect(() => {
        itemsContainerApi.start(
            {
                height: `${open ? itemsContainerRef.current?.offsetHeight : 0}px`,
                opacity: open ? 1 : 0,
                config: {
                    easing: theme.transitions.easing.quintFn,
                    duration: theme.transitions.duration.quint,
                },
            }
        )
    }, [open, isWindowResizing])

    return (
        <Box
            sx={{

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
                    className='drawer-section'
                    variant='pMediumBoldWt'
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

            <animated.div
                style={{
                    overflow: 'hidden',
                    ...itemsContainerStyle,
                }}
            >
                <Box
                    ref={itemsContainerRef}
                >
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
            </animated.div >
        </Box >
    )
}