import * as React from 'react';

// TODO
// https://stackoverflow.com/questions/16302483/event-to-detect-when-positionsticky-is-triggered
export default function useIsStuck(
    ref: React.RefObject<HTMLElement>,
    observerSettings: IntersectionObserverInit = {
        rootMargin: '-1px 0px 0px 0px',
        root: null,
        threshold: [1],
    }
) {
    const [isStuck, setIsStuck] = React.useState(false)

    React.useEffect(() => {
        const refEl = ref.current
        const observer = new IntersectionObserver(
            ([e]) => setIsStuck(e.intersectionRatio < 1),
            observerSettings,
        )

        if (refEl == null) {
            // TODO: log something
            return
        }
        observer.observe(refEl)

        return () => {
            observer.unobserve(refEl)
        }
    }, [])

    return isStuck
}