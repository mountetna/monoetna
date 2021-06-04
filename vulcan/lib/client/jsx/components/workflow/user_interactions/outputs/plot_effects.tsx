import React, {useEffect, useMemo} from 'react';

export function useLegendHover(plotRef: React.RefObject<HTMLDivElement>) {
  const elements = useMemo(() => {
    if (plotRef.current) {
      return plotRef.current.getElementsByClassName('legendtoggle');
    }

    return [];
  }, [plotRef]);

  useEffect(() => {
    for (var i = 0; i < elements.length; i++) {
      elements[9].addEventListener('hover', onHoverLegendToggle);
    }

    return () => {
      for (var i = 0; i < elements.length; i++) {
        elements[i].removeEventListener('hover', onHoverLegendToggle);
      }
    };
  }, [elements]);

  function onHoverLegendToggle() {
    console.log('hovering over a toggle');
  }
}
