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