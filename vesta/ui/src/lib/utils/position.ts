import * as React from 'react'
import _ from 'lodash'


interface TwoDCoords {
    x: number
    y: number
}

export function usePointerPosition(mode: 'pixels' | 'fraction' = 'fraction', throttleIntervalMs: number = 100): TwoDCoords | null {
    const [coords, setCoords] = React.useState<TwoDCoords | null>(null)

    let handlePointerMove = (event: PointerEvent) => {
        const viewportX = event.clientX
        const viewportY = event.clientY
        const viewportWidth = window.innerWidth
        const viewportHeight = window.innerHeight

        setCoords({
            x: mode === 'fraction' ? viewportX / viewportWidth : viewportX,
            y: mode === 'fraction' ? viewportY / viewportHeight : viewportY,
        })
    }

    handlePointerMove = _.throttle(
        handlePointerMove,
        throttleIntervalMs,
        {
            leading: true,
            trailing: true,
        },
    )

    React.useEffect(() => {
        document.addEventListener('pointermove', handlePointerMove)
        return () => document.removeEventListener('pointermove', handlePointerMove)
    }, [])

    return coords
}

interface ThreeDCoords {
    x: number
    y: number
    z: number
}

export function useDeviceOrientation(throttleIntervalMs: number = 100): ThreeDCoords | null {
    const [coords, setCoords] = React.useState<ThreeDCoords | null>(null)

    let handleDeviceOrientation = (event: DeviceOrientationEvent) => {
        if (event.alpha == null || event.beta == null || event.gamma == null) {
            return
        }

        setCoords({
            x: -event.beta,
            y: event.gamma,
            z: -event.alpha,
        })
    }

    handleDeviceOrientation = _.throttle(
        handleDeviceOrientation,
        throttleIntervalMs,
        {
            leading: true,
            trailing: true,
        },
    )

    React.useEffect(() => {
        window.addEventListener('deviceorientation', handleDeviceOrientation)
        return () => window.removeEventListener('deviceorientation', handleDeviceOrientation)
    }, [])

    return coords
}