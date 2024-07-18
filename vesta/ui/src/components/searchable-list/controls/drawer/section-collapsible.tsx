'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Typography from '@mui/material/Typography';
import ButtonBase from '@mui/material/ButtonBase';
import { useTheme } from '@mui/material';
import { useSpring, animated } from '@react-spring/web';
import _ from 'lodash'

import { useWindowDimensions } from '@/lib/utils/responsive';
import { DrawerSectionProps } from './models';


const activeClass = 'active'

export default function DrawerSectionCollapsible({
    name,
    items,
    activeKeys,
    onClickItem,
}: DrawerSectionProps) {
    const theme = useTheme()

    const [open, setOpen] = React.useState(false)
    const { isResizing: isWindowResizing } = useWindowDimensions()

    // Manage open/close
    const rootRef = React.useRef<HTMLElement>()
    const [rootStyle, rootApi] = useSpring(() => ({
        height: '0px',
        opacity: 0,
    }), [open])

    React.useEffect(() => {
        rootApi.start(
            {
                height: `${open ? rootRef.current?.offsetHeight : 0}px`,
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
            <Typography
                className='drawer-section'
                variant='pMediumBoldWt'
            >
                {_.startCase(name)}
            </Typography>

            {items.map(item => (
                <Typography
                    key={item.key}
                    className={`drawer-item${activeKeys.has(item.key) ? ` ${activeClass}` : ''}`}
                    variant='pMedium'
                >
                    <ButtonBase
                        onClick={() => onClickItem(item)}
                        sx={{

                        }}
                    >
                        {item.label}
                    </ButtonBase>
                </Typography>
            ))}
        </Box>
    )
}