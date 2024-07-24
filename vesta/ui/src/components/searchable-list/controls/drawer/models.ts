export interface DrawerItem {
    label: string
    key: string
    type: string
}

export type DisplayStyle = 'default' | 'expandable'

export interface DrawerSectionProps {
    name: string;
    items: DrawerItem[];
    activeKeys: Set<string>;
    onClickItem: (item: DrawerItem) => void;
}

export enum DrawerMainContentClass {
    base = 'drawer-main-content',
    default = 'drawer-main-content-default',
    expandable = 'drawer-main-content-expandable',
}

export enum DrawerSectionClass {
    base = 'drawer-section',
    viewSets = 'drawer-section-view-sets',
    default = 'drawer-section-default',
    expandable = 'drawer-section-expandable',
}