import * as React from 'react'
import _ from 'lodash'


interface Coords {
    x: number
    y: number
}

export function useMousePosition(mode: 'pixels' | 'fraction' = 'fraction', throttleIntervalMs: number = 100): Coords {
    const [coords, setCoords] = React.useState<Coords>({ x: 0, y: 0 })

    let handleMouseMove = (event: MouseEvent) => {
        const viewportX = event.clientX
        const viewportY = event.clientY
        const viewportWidth = window.innerWidth
        const viewportHeight = window.innerHeight

        setCoords({
            x: mode === 'fraction' ? viewportX / viewportWidth : viewportX,
            y: mode === 'fraction' ? viewportY / viewportHeight : viewportY,
        })
    }

    handleMouseMove = _.throttle(
        handleMouseMove,
        throttleIntervalMs,
        {
            leading: true,
            trailing: true,
        },
    )

    React.useEffect(() => {
        document.addEventListener('mousemove', handleMouseMove)
        return () => document.removeEventListener('mousemove', handleMouseMove)
    }, [])

    return coords
}