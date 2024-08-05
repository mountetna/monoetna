'use client'

import * as React from 'react'
import Box from '@mui/system/Box'
import Collapse from '@mui/material/Collapse'
import Fade from '@mui/material/Fade'
import { useTheme } from '@mui/material';

import { DrawerItem, DisplayStyle, DrawerSectionProps, DrawerSectionClass, DrawerMainContentClass } from './models';
import DrawerSectionDefault from './section-default';
import DrawerSectionExpandable from './section-expandable';
import Button from '@/components/inputs/button';
import { TransitionProps } from '@mui/material/transitions'


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
    const [sectionOpens, setSectionOpens] = React.useState(Object.keys(drawerItemsByType).map(_ => false).concat(false))

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
                break
            case 'expandable':
                drawerSection = (
                    <DrawerSectionExpandable
                        key={sectionName}
                        open={sectionOpens[sectionOpenIdx]}
                        onSetOpen={(open) => handleSetSectionOpen(open, sectionOpenIdx)}
                        {...sectionProps}
                    />
                )
                break
        }

        viewSets = (
            <Box
                key={sectionName}
                className={[DrawerSectionClass.base, DrawerSectionClass.viewSets].join(' ')}
            >
                {drawerSection}
            </Box>
        )
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
                    // ref={rootRef}
                    className='drawer-main-content-container'
                >
                    <Box
                        className={[DrawerMainContentClass.base, displayStyle === 'default' ? DrawerMainContentClass.default : DrawerMainContentClass.expandable].join(' ')}
                        sx={{
                            display: 'flex',
                            borderRadius: '30px',
                            bgcolor: 'utilityWhite.main',
                            [`& .${DrawerSectionClass.base}`]: {
                            },
                            [`&.${DrawerMainContentClass.default}`]: {
                                flexDirection: 'row',
                                gap: '16px',
                                p: '24px',
                                minWidth: 'fit-content',
                                '& > *:not(:last-child)': {
                                    borderRight: `1px solid ${theme.palette.ground.grade75}`,
                                },
                                [`& .${DrawerSectionClass.viewSets}`]: {
                                },
                            },
                            [`& .${DrawerSectionClass.default}`]: {
                            },
                            [`&.${DrawerMainContentClass.expandable}`]: {
                                flexDirection: 'column',
                                gap: '16px',
                                p: '24px',
                                [`& .${DrawerSectionClass.viewSets}`]: {
                                    pb: '10px',
                                },
                            },
                            [`& .${DrawerSectionClass.expandable}`]: {
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
                            let drawerSectionClass
                            switch (displayStyle) {
                                case 'default':
                                    drawerSection = (
                                        <DrawerSectionDefault
                                            {...sectionProps}
                                        />
                                    )
                                    drawerSectionClass = DrawerSectionClass.default
                                    break
                                case 'expandable':
                                    drawerSection = (
                                        <DrawerSectionExpandable
                                            open={sectionOpens[index]}
                                            onSetOpen={(open) => handleSetSectionOpen(open, index)}
                                            {...sectionProps}
                                        />
                                    )
                                    drawerSectionClass = DrawerSectionClass.expandable
                                    break
                            }

                            return (
                                <Box
                                    key={section}
                                    className={[DrawerSectionClass.base, drawerSectionClass].join(' ')}
                                >
                                    {drawerSection}
                                </Box>
                            )
                        })}

                        {showButton && <Button
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
                        />}
                    </Box>
                </Box>
            </Fade>
        </Collapse>
    )
}