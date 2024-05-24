import * as React from 'react';

// This works for elements not already at the top of page.
// Probably not a generalized solution.
export default function useIsStuck(
    ref: React.RefObject<HTMLElement>,
    observerSettings: IntersectionObserverInit = {
        root: null,
        threshold: [1],
    }
) {
    const [isStuck, setIsStuck] = React.useState(false)

    React.useEffect(() => {
        const refEl = ref.current
        
        if (document == undefined || refEl == null) {
            return
        }

        const detectorEl = document.createElement('div')
        detectorEl.style.cssText = 'visibility: hidden'

        document.body.insertBefore(detectorEl, refEl)
        
        const observer = new IntersectionObserver(
            ([e]) => setIsStuck(e.intersectionRatio < 1),
            observerSettings,
        )
        
        observer.observe(detectorEl)

        return () => {
            observer.unobserve(detectorEl)
            detectorEl.remove()
        }
    }, [])

    return isStuck
}