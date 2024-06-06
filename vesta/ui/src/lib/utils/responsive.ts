import { useState, useEffect } from 'react'
import _ from 'lodash'


export function useWindowDimensions(triggerDelayMs: number = 100) {
    const [dimensions, setDimensions] = useState<number[]>()
    const [isResizing, setIsResizing] = useState<boolean>(false)

    useEffect(() => {
        const dimensionsHandler = _.throttle(() => {
            setDimensions([window.innerWidth, window.innerHeight])
        }, triggerDelayMs)

        const resizingTrueHandler = _.throttle(() => {
            setIsResizing(true)
        }, triggerDelayMs)

        const resizingFalseHandler = _.debounce(() => {
            resizingTrueHandler.cancel()
            setIsResizing(false)
        }, triggerDelayMs)

        const handlers = [
            dimensionsHandler,
            resizingTrueHandler,
            resizingFalseHandler,
        ]

        for (const handler of handlers) {
            window.addEventListener('resize', handler)
        }

        return () => {
            for (const handler of handlers) {
                window.removeEventListener('resize', handler)
            }
        }
    }, [triggerDelayMs])

    return { dimensions, isResizing }
}