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

export enum DrawerClasses {
    root = 'drawer-main-content-container',
    mainContent = 'drawer-main-content',
    extraControls = 'drawer-extra-controls-container',
}

export enum DrawerSectionsContainerClasses {
    base = 'drawer-sections-container',
    default = 'drawer-sections-container-default',
    expandable = 'drawer-sections-container-expandable',
}

export enum DrawerSectionClasses {
    base = 'drawer-section',
    viewSets = 'drawer-section-view-sets',
    default = 'drawer-section-default',
    expandable = 'drawer-section-expandable',
}