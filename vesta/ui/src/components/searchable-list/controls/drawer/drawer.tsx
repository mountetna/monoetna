'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import { useTheme } from '@mui/material';
import { useSpring, animated, useIsomorphicLayoutEffect } from '@react-spring/web';

import { useWindowDimensions } from '@/lib/utils/responsive';
import { DrawerItem, DisplayStyle, DrawerSectionProps } from './models';
import DrawerSectionDefault from './section-default';
import DrawerSectionExpandable from './section-expandable';


// TODO: make item a generic and
// add getters like getItemKey
export default function Drawer<Item>({
    items,
    getItemLabel,
    getItemKey,
    getItemType,
    activeItems,
    onChange,
    open,
    displayStyle = 'default',
}: {
    items: Item[],
    getItemLabel: (item: Item) => string,
    getItemKey: (item: Item) => string,
    getItemType: (item: Item) => string,
    activeItems: Item[],
    onChange: (activeItems: Item[]) => void,
    open: boolean,
    displayStyle?: DisplayStyle,
}) {
    const theme = useTheme()
    const [finishedMountAnimation, setFinishedMountAnimation] = React.useState(false)

    const itemsByKey: Record<string, Item> = {}
    const drawerItemsByKey: Record<string, DrawerItem> = {}
    const drawerItemsByType: Record<string, DrawerItem[]> = {}

    items.forEach(item => {
        const drawerItem = {
            label: getItemLabel(item),
            key: getItemKey(item),
            type: getItemType(item),
        } as DrawerItem

        itemsByKey[drawerItem.key] = item
        drawerItemsByKey[drawerItem.key] = drawerItem

        if (!(drawerItem.type in drawerItemsByType)) {
            drawerItemsByType[drawerItem.type] = [drawerItem]
            return
        }
        drawerItemsByType[drawerItem.type].push(drawerItem)
    })

    const activeKeys = new Set(activeItems.map(item => getItemKey(item)))

    // Manage open/close
    const { isResizing: isWindowResizing } = useWindowDimensions()

    const [sectionOpens, setSectionOpens] = React.useState(Object.keys(drawerItemsByType).map(_ => false))

    const handleSetSectionOpen = (open: boolean, sectionIndex: number) => {
        const newSectionOpens = [...sectionOpens]
        newSectionOpens[sectionIndex] = open
        setSectionOpens(newSectionOpens)
    }

    const rootRef = React.useRef<HTMLElement>()
    const [rootStyle, rootApi] = useSpring(() => ({
        height: '0px',
        opacity: 0,
    }), [open])

    useIsomorphicLayoutEffect(() => {
        if (!open && finishedMountAnimation) {
            rootApi.set({ height: `${(rootRef.current?.offsetHeight || 0)}px`, })
            rootApi.set({ height: `${(rootRef.current?.offsetHeight || 0)}px`, })
        }
        rootApi.start(
            {
                height: `${open ? (rootRef.current?.offsetHeight || 0) : 0}px`,
                opacity: open ? 1 : 0,
                config: {
                    easing: theme.transitions.easing.quintFn,
                    duration: theme.transitions.duration.quint,
                },
                onStart: () => {
                    const animatedEl = rootRef.current?.parentElement
                    if (!open && animatedEl) {
                        animatedEl.className = ''
                    }
                },
                onRest: () => {
                    const animatedEl = rootRef.current?.parentElement
                    if (open && animatedEl) {
                        animatedEl.className = 'full-height'
                    }

                    setFinishedMountAnimation(true)
                }
            }
        )
    }, [open, isWindowResizing])

    const handleClickItem = (item: DrawerItem) => {
        let newActiveItems: Item[]

        if (activeKeys.has(item.key)) {
            newActiveItems = activeItems
                .filter(_item => getItemKey(_item) !== item.key)
        } else {
            newActiveItems = [...activeItems, itemsByKey[item.key]]
        }

        return onChange(newActiveItems)
    }

    return (
        <Box
            sx={{
                '& > *.full-height': {
                    height: 'auto !important',
                },
            }}
        >
            <animated.div
                style={{
                    overflow: 'hidden',
                    ...rootStyle,
                }}
            >
                <Box
                    ref={rootRef}
                    sx={{
                        p: '24px',
                        borderRadius: '30px',
                        bgcolor: 'utilityWhite.main',
                    }}
                >
                    {/* TODO: make separate container components? */}
                    {Object.entries(drawerItemsByType).map(([section, items], index) => {
                        const sectionProps: DrawerSectionProps = {
                            name: section,
                            items: items,
                            activeKeys: activeKeys,
                            onClickItem: (item) => handleClickItem(item),
                        }

                        switch (displayStyle) {
                            case 'default':
                                return (
                                    <DrawerSectionDefault
                                        key={section}
                                        {...sectionProps}
                                    />
                                )
                            case 'expandable':
                                return (
                                    <DrawerSectionExpandable
                                        key={section}
                                        open={sectionOpens[index]}
                                        onSetOpen={(open) => handleSetSectionOpen(open, index)}
                                        {...sectionProps}
                                    />
                                )
                        }
                    })}

                    <Box>
                        Export placeholder
                    </Box>
                </Box>
            </animated.div>
        </Box>
    )
}