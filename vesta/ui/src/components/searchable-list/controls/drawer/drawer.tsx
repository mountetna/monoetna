'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import { useTheme } from '@mui/material';
import { useSpring, animated } from '@react-spring/web';

import { useWindowDimensions } from '@/lib/utils/responsive';
import { DrawerItem, DisplayStyle, DrawerSectionProps } from './models';
import DrawerSectionDefault from './section-default';
import DrawerSectionCollapsible from './section-collapsible';


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

    const itemsByKey: Record<string, Item> = {}
    const drawerItemsByKey: Record<string, DrawerItem> = {}
    items.forEach(item => {
        const drawerItem = {
            label: getItemLabel(item),
            key: getItemKey(item),
            type: getItemType(item),
        } as DrawerItem

        itemsByKey[drawerItem.key] = item
        drawerItemsByKey[drawerItem.key] = drawerItem

        return drawerItem
    })

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

    const itemsByType: Record<string, DrawerItem[]> = {}
    for (const item of Object.values(drawerItemsByKey)) {
        if (!(item.type in itemsByType)) {
            itemsByType[item.type] = [item]
            continue
        }

        itemsByType[item.type].push(item)
    }

    const activeKeys = new Set(activeItems.map(item => getItemKey(item)))

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
        <animated.div
            style={{
                overflow: 'hidden',
                ...rootStyle
            }}
        >
            <Box
                ref={rootRef}
            >
                {Object.entries(itemsByType).map(([section, items]) => {
                    const sectionProps: DrawerSectionProps = {
                        name: section,
                        items: items,
                        activeKeys: activeKeys,
                        onClickItem: (item) => handleClickItem(item),
                    }

                    switch (displayStyle) {
                        case 'default':
                            return <DrawerSectionDefault key={section} {...sectionProps} />
                        case 'collapsible':
                            return <DrawerSectionCollapsible key={section} {...sectionProps} />
                    }
                })}
            </Box>
        </animated.div>
    )
}