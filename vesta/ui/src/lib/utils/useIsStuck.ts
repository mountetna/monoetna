import * as React from 'react';

// TODO
// https://stackoverflow.com/questions/16302483/event-to-detect-when-positionsticky-is-triggered
export default function useIsStuck(el: HTMLElement) {
    const [isStuck, setIsStuck] = React.useState(false)

    React.useEffect(() => {

    }, [])

    return isStuck
}