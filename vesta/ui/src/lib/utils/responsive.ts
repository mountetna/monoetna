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

        window.addEventListener('resize', dimensionsHandler)
        window.addEventListener('resize', resizingTrueHandler)
        window.addEventListener('resize', resizingFalseHandler)

        return () => {
            window.removeEventListener('resize', dimensionsHandler)
            window.removeEventListener('resize', resizingTrueHandler)
            window.removeEventListener('resize', resizingFalseHandler)
        }
    }, [triggerDelayMs])

    return { dimensions, isResizing }
}