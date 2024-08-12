'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Collapse from '@mui/material/Collapse'
import Fade from '@mui/material/Fade'
import { useTheme } from '@mui/material';

import { DrawerItem, DisplayStyle, DrawerSectionProps, DrawerSectionClasses, DrawerSectionsContainerClasses, DrawerClasses } from './models';
import DrawerSectionDefault from './section-default';
import DrawerSectionExpandable from './section-expandable';
import Button from '@/components/inputs/button';
import { TransitionProps } from '@mui/material/transitions'
import ToggleGroup, { Classes as ToggleGroupClasses } from '@/components/inputs/toggle-group'
import { useBreakpoint } from '@/lib/utils/responsive'


export default function Drawer<Item>({
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

    const breakpoint = useBreakpoint()

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
    // Section open for each item type + View Set section
    const [sectionOpens, setSectionOpens] = React.useState(Object.keys(drawerItemsByType).map(_ => false).concat([false]))

    const handleSetSectionOpen = (open: boolean, sectionIndex: number) => {
        const newSectionOpens = [...sectionOpens]
        newSectionOpens[sectionIndex] = open
        setSectionOpens(newSectionOpens)
    }

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
    let mobileViewSets
    let tabletDesktopViewSets
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
                break
            case 'expandable':
                drawerSection = (
                    <DrawerSectionExpandable
                        key={sectionName}
                        open={sectionOpens[sectionOpenIdx]}
                        onSetOpen={(open) => handleSetSectionOpen(open, sectionOpenIdx)}
                        variant='withBackground'
                        {...sectionProps}
                    />
                )
                break
        }

        mobileViewSets = (
            <Box
                key={sectionName}
                className={[DrawerSectionClasses.base, DrawerSectionClasses.viewSets].join(' ')}
                sx={{
                    [theme.breakpoints.up('tablet')]: {
                        display: 'none',
                    },
                }}
            >
                {drawerSection}
            </Box>
        )

        let activeViewSetItemIndex = 0
        const activeViewSetItemKey = getItemKey(activeViewSetItem)
        for (const [idx, item] of viewSetDrawerItems.entries()) {
            if (activeViewSetItemKey === item.key) {
                activeViewSetItemIndex = idx
                break
            }
        }

        tabletDesktopViewSets = (
            <Box
                sx={{
                    [`& .${ToggleGroupClasses.root}`]: {
                        bgcolor: 'utilityHighlight.main',
                        borderRadius: '40px',
                        [`& .${ToggleGroupClasses.toggleRoot}:not(.${ToggleGroupClasses.toggleSelected})`]: {
                            bgcolor: 'utilityWhite.main',
                        },
                    },
                }}
            >
                <ToggleGroup
                    valueIdx={activeViewSetItemIndex}
                    values={viewSetDrawerItems.map(val => val.label)}
                    onChange={(idx) => onChangeViewSet(viewSetItems[idx])}
                />
            </Box>
        )
    }

    let button
    let mobileButton
    let tabletDesktopButton
    if (showButton) {
        button = (
            <Button
                label={buttonLabel}
                onClick={onClickButton}
                strokeVariant='stroked'
                sizeVariant='large'
                typographyVariant='pBodyMediumWt'
                sx={{
                    width: 'auto',
                    height: 'auto',
                }}
                disabled={buttonDisabled}
            />
        )

        if (breakpoint === 'mobile') {
            mobileButton = button
        } else {
            tabletDesktopButton = button
        }
    }

    const animationProps: TransitionProps = {
        in: open,
        easing: theme.transitions.easing.quint,
        timeout: theme.transitions.duration.quint,
    }

    return (
        <Collapse
            {...animationProps}
        >
            <Fade
                {...animationProps}
            >
                <Box
                    className={DrawerClasses.root}
                >
                    <Box
                        className={DrawerClasses.mainContent}
                        sx={{
                            borderRadius: '30px',
                            bgcolor: 'utilityWhite.main',
                            [theme.breakpoints.up('tablet')]: {
                                display: 'flex',
                                flexDirection: 'column',
                                gap: '24px',
                                p: '16px',
                                pb: '24px',
                            },
                        }}
                    >
                        <Box
                            className={DrawerClasses.extraControls}
                            sx={{
                                display: 'none',
                                [theme.breakpoints.up('tablet')]: {
                                    display: 'flex',
                                    flexDirection: 'row',
                                    justifyContent: 'space-between',
                                    alignItems: 'center',
                                },
                            }}
                        >
                            {tabletDesktopViewSets}

                            {tabletDesktopButton}
                        </Box>

                        <Box
                            className={[DrawerSectionsContainerClasses.base, displayStyle === 'default' ? DrawerSectionsContainerClasses.default : DrawerSectionsContainerClasses.expandable].join(' ')}
                            sx={{
                                display: 'flex',
                                [`&.${DrawerSectionsContainerClasses.base}`]: {
                                    gap: '16px',
                                    p: '24px',
                                    [theme.breakpoints.up('tablet')]: {
                                        p: '0px',
                                    },
                                },
                                [`& .${DrawerSectionClasses.base}`]: {
                                },
                                [`&.${DrawerSectionsContainerClasses.default}`]: {
                                    flexDirection: 'row',
                                    minWidth: 'fit-content',
                                    '& > *:not(:last-child)': {
                                        borderRight: `1px solid ${theme.palette.ground.grade75}`,
                                    },
                                    [`& .${DrawerSectionClasses.viewSets}`]: {
                                    },
                                },
                                [`& .${DrawerSectionClasses.default}`]: {
                                },
                                [`&.${DrawerSectionsContainerClasses.expandable}`]: {
                                    flexDirection: 'column',
                                    [`& .${DrawerSectionClasses.viewSets}`]: {
                                    },
                                },
                                [`& .${DrawerSectionClasses.expandable}`]: {
                                    '&:last-of-type': {
                                        mb: '8px',
                                    },
                                },
                            }}
                        >
                            {mobileViewSets}

                            {/* TODO: make separate container components? */}
                            {Object.entries(drawerItemsByType).map(([section, items], index) => {
                                const sectionProps: DrawerSectionProps = {
                                    name: section,
                                    items: items,
                                    activeKeys,
                                    onClickItem: (item) => handleClickItem(item),
                                }

                                let drawerSection
                                let drawerSectionClass
                                switch (displayStyle) {
                                    case 'default':
                                        drawerSection = (
                                            <DrawerSectionDefault
                                                {...sectionProps}
                                            />
                                        )
                                        drawerSectionClass = DrawerSectionClasses.default
                                        break
                                    case 'expandable':
                                        drawerSection = (
                                            <DrawerSectionExpandable
                                                open={sectionOpens[index]}
                                                onSetOpen={(open) => handleSetSectionOpen(open, index)}
                                                variant='default'
                                                {...sectionProps}
                                            />
                                        )
                                        drawerSectionClass = DrawerSectionClasses.expandable
                                        break
                                }

                                return (
                                    <Box
                                        key={section}
                                        className={[DrawerSectionClasses.base, drawerSectionClass].join(' ')}
                                    >
                                        {drawerSection}
                                    </Box>
                                )
                            })}

                            {mobileButton}
                        </Box>
                    </Box>
                </Box>
            </Fade>
        </Collapse>
    )
}