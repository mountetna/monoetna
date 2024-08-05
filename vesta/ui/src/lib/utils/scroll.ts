import { getMainNavHeight } from '@/components/nav/main-nav';
import { Breakpoint } from '@mui/material';


export interface ScrollToProps {
    top?: number
    left?: number
    ignoreMainNav?: boolean
    behavior?: ScrollBehavior
}

export function scrollTo(props: ScrollToProps, breakpoint: Breakpoint) {
    const ignoreMainNav = props.ignoreMainNav === undefined ? false : props.ignoreMainNav
    const topCompensation = ignoreMainNav ? 0 : getMainNavHeight(breakpoint)

    window.scrollTo({
        top: props.top !== undefined ? props.top - topCompensation : props.top,
        left: props.left,
        behavior: props.behavior === undefined ? 'smooth' : props.behavior,
    })
}