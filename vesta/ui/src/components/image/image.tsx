import * as React from 'react'
import BaseImage, { ImageProps } from 'next/image'
import { useTheme } from '@mui/material'


interface Props extends ImageProps {
    hideBeforeLoad?: boolean
}


export default function Image(props: Props) {
    const {
        hideBeforeLoad,
        ...imageProps
    } = props

    const theme = useTheme()

    const [loaded, setLoaded] = React.useState(hideBeforeLoad === undefined ? true : hideBeforeLoad)

    return (
        <BaseImage
            {...imageProps}
            onLoad={event => {
                setLoaded(true)
                props.onLoad && props.onLoad(event)
            }}
            style={{
                ...(props.style ? props.style : {}),
                opacity: loaded ? 1 : 0,
                transition: theme.transitions.create(
                    ['opacity'],
                    {
                        easing: theme.transitions.easing.swell,
                        duration: theme.transitions.duration.swell,
                    },
                ),
            }}
        />
    )
}