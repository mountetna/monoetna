import * as React from 'react';

// This works for elements not already at the top of page.
// Probably not a generalized solution.
// TODO: fix glitches https://css-tricks.com/how-to-create-a-shrinking-header-on-scroll-without-javascript/
// https://stackoverflow.com/questions/62249752/why-does-changing-the-height-of-a-positionsticky-element-alter-the-scroll-pos
export default function useIsStuck(
    ref: React.RefObject<HTMLElement>,
    observerSettings: IntersectionObserverInit = {
        root: null,
        threshold: [0],
    }
) {
    const [isStuck, setIsStuck] = React.useState(false)

    React.useEffect(() => {
        const refEl = ref.current
        const parentEl = ref.current?.parentNode

        if (document == undefined || refEl == null || parentEl == null) {
            return
        }

        const detectorEl = document.createElement('div')
        detectorEl.style.visibility = 'hidden'

        parentEl.insertBefore(detectorEl, refEl)

        const observer = new IntersectionObserver(
            ([e]) => setIsStuck(!e.isIntersecting),
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