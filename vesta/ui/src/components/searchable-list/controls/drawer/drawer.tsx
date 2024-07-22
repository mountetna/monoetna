'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import { useTheme } from '@mui/material';
import { useSpring, animated, useIsomorphicLayoutEffect } from '@react-spring/web';

import { useWindowDimensions } from '@/lib/utils/responsive';
import { DrawerItem, DisplayStyle, DrawerSectionProps } from './models';
import DrawerSectionDefault from './section-default';
import DrawerSectionExpandable from './section-expandable';
import Button from '@/components/inputs/button';


const baseDrawerSectionClass = 'drawer-section'
const viewSetsSectionClass = 'drawer-section-view-sets'
const defaultSectionClass = 'drawer-section-default'


// TODO: make item a generic and
// add getters like getItemKey
export default function Drawer<Item, ViewSetItem>({
    items,
    viewSetItems,
    getItemLabel,
    getItemKey,
    getItemType,
    activeItems,
    onChange,
    activeViewSetItem,
    onChangeViewSet,
    showViewSets,
    showButton,
    buttonLabel,
    onClickButton,
    buttonDisabled,
    open,
    displayStyle = 'default',
}: {
    items: Item[],
    viewSetItems: Item[],
    getItemLabel: (item: Item) => string,
    getItemKey: (item: Item) => string,
    getItemType: (item: Item) => string,
    activeItems: Item[],
    onChange: (activeItems: Item[]) => void,
    activeViewSetItem: Item,
    onChangeViewSet: (activeViewSetItem: Item) => void,
    showViewSets: boolean,
    showButton: boolean,
    buttonLabel: string,
    onClickButton: () => void,
    buttonDisabled?: boolean,
    open: boolean,
    displayStyle?: DisplayStyle,
}) {
    const theme = useTheme()
    const [finishedMountAnimation, setFinishedMountAnimation] = React.useState(false)

    const itemsByKey: Record<string, Item> = {}
    const drawerItemsByType: Record<string, DrawerItem[]> = {}
    const viewSetItemsByKey: Record<string, Item> = {}
    const viewSetDrawerItems: DrawerItem[] = []

    items.forEach(item => {
        const drawerItem: DrawerItem = {
            label: getItemLabel(item),
            key: getItemKey(item),
            type: getItemType(item),
        }

        itemsByKey[drawerItem.key] = item

        if (!(drawerItem.type in drawerItemsByType)) {
            drawerItemsByType[drawerItem.type] = [drawerItem]
            return
        }
        drawerItemsByType[drawerItem.type].push(drawerItem)
    })
    viewSetItems.forEach(item => {
        const drawerItem: DrawerItem = {
            label: getItemLabel(item),
            key: getItemKey(item),
            type: getItemType(item),
        }

        viewSetItemsByKey[drawerItem.key] = item
        viewSetDrawerItems.push(drawerItem)
    })

    const activeKeys = new Set(activeItems.map(item => getItemKey(item)))

    // Manage open/close
    const { isResizing: isWindowResizing } = useWindowDimensions()

    // Section open for each item type + View Set section
    const [sectionOpens, setSectionOpens] = React.useState(Object.keys(drawerItemsByType).map(_ => false).concat(false))

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

    const handleClickViewSetItem = (item: DrawerItem) => {
        onChangeViewSet(viewSetItemsByKey[item.key])
    }

    // Manage view sets
    let viewSets
    if (showViewSets) {
        const sectionName = 'view'
        const sectionOpenIdx = sectionOpens.length - 1

        const sectionProps: DrawerSectionProps = {
            name: sectionName,
            items: viewSetDrawerItems,
            activeKeys: new Set([getItemKey(activeViewSetItem)]),
            onClickItem: handleClickViewSetItem
        }

        let drawerSection
        switch (displayStyle) {
            case 'default':
                drawerSection = (
                    <DrawerSectionDefault
                        key={sectionName}
                        {...sectionProps}
                    />
                )
            case 'expandable':
                drawerSection = (
                    <DrawerSectionExpandable
                        key={sectionName}
                        open={sectionOpens[sectionOpenIdx]}
                        onSetOpen={(open) => handleSetSectionOpen(open, sectionOpenIdx)}
                        {...sectionProps}
                    />
                )
        }

        viewSets = (
            <Box
                key={sectionName}
                className={[baseDrawerSectionClass, viewSetsSectionClass].join(' ')}
            >
                {drawerSection}
            </Box>
        )
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
                        display: 'flex',
                        flexDirection: 'column',
                        gap: '16px',
                        p: '24px',
                        borderRadius: '30px',
                        bgcolor: 'utilityWhite.main',
                        [`& .${baseDrawerSectionClass}`]: {
                        },
                        [`& .${viewSetsSectionClass}`]: {
                            pb: '10px',
                        },
                        [`& .${defaultSectionClass}`]: {
                            pb: '16px',
                            borderBottom: `1px solid ${theme.palette.ground.grade75}`,
                            '&:last-of-type': {
                                mb: '8px',
                            },
                        },
                    }}
                >
                    {viewSets}

                    {/* TODO: make separate container components? */}
                    {Object.entries(drawerItemsByType).map(([section, items], index) => {
                        const sectionProps: DrawerSectionProps = {
                            name: section,
                            items: items,
                            activeKeys,
                            onClickItem: (item) => handleClickItem(item),
                        }

                        let drawerSection
                        switch (displayStyle) {
                            case 'default':
                                drawerSection = (
                                    <DrawerSectionDefault
                                        {...sectionProps}
                                    />
                                )
                            case 'expandable':
                                drawerSection = (
                                    <DrawerSectionExpandable
                                        open={sectionOpens[index]}
                                        onSetOpen={(open) => handleSetSectionOpen(open, index)}
                                        {...sectionProps}
                                    />
                                )
                        }

                        return (
                            <Box
                                key={section}
                                className={[baseDrawerSectionClass, defaultSectionClass].join(' ')}
                            >
                                {drawerSection}
                            </Box>
                        )
                    })}

                    {showButton &&
                        <Button
                            label={buttonLabel}
                            onClick={onClickButton}
                            variant='large'
                            typographyVariant='pBodyMediumWt'
                            sx={{
                                width: 'auto',
                                height: 'auto',
                                px: '16px',
                                py: '8px',
                                borderRadius: '10px',
                            }}
                            disabled={buttonDisabled}
                        />
                    }
                </Box>
            </animated.div>
        </Box>
    )
}