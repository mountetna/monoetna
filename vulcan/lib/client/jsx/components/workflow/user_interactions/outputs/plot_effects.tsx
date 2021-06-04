import React, {useCallback, useEffect, useMemo} from 'react';

export function useLegendHover(plot: HTMLDivElement | null) {
  const elements = useMemo(() => {
    if (plot) {
      return plot.getElementsByClassName('legendtoggle');
    }

    return [];
  }, [plot]);

  const dimOtherTraces = useCallback(
    (element) => {
      // Set opacity of the other traces to 0.5
      for (var i = 0; i < elements.length; i++) {
        if (elements[i] === element) continue;

        console.log('dimming i', i);
      }
    },
    [elements]
  );

  const onHoverLegendToggle = useCallback(
    (e: any) => {
      console.log(e);
      console.log('hovering over a toggle');
      dimOtherTraces(e.target);
    },
    [dimOtherTraces]
  );

  useEffect(() => {
    if (elements) {
      for (var i = 0; i < elements.length; i++) {
        elements[i].addEventListener('mouseover', onHoverLegendToggle);
      }
    }

    return () => {
      for (var i = 0; i < elements.length; i++) {
        elements[i].removeEventListener('mouseover', onHoverLegendToggle);
      }
    };
  }, [elements, onHoverLegendToggle]);

  if (!plot) return;
}
