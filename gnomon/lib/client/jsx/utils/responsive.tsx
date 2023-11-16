import React, { useState, useEffect } from 'react';
import _ from 'lodash';


export function useWindowDimensions(triggerDelayMs: number = 100) {
    const [dimensions, setDimensions] = useState([
        window.innerWidth,
        window.innerHeight,
    ]);

    useEffect(() => {
        const debouncedResizeHandler = _.debounce(() => {
            setDimensions([window.innerWidth, window.innerHeight]);
        }, triggerDelayMs);

        window.addEventListener('resize', debouncedResizeHandler);

        return () => window.removeEventListener('resize', debouncedResizeHandler);
    }, []);

    return dimensions;
}