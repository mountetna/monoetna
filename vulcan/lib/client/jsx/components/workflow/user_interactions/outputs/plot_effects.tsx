import React, {useCallback, useEffect, useMemo} from 'react';

export function useLegendHover(
  plot: HTMLDivElement | null,
  onHover: (index: number | undefined) => void
) {
  const elements = useMemo(() => {
    if (plot) {
      return plot.getElementsByClassName('legendtoggle');
    }

    return [];
  }, [plot]);

  const onHoverLegendToggle = useCallback(
    (e: any) => {
      let index;

      for (var i = 0; i < elements.length; i++) {
        if (elements[i] !== e.target) continue;

        index = i;
        break;
      }
      onHover(index);
    },
    [onHover, elements]
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
