import React, {useCallback, useMemo, useState} from 'react';

export function useLegendHover(
  plot: HTMLDivElement | null,
  onHover: React.Dispatch<React.SetStateAction<number>>
): {attachEventListeners: () => void; removeEventListeners: () => void} {
  const [listenersAttached, setListenersAttached] = useState(false);

  const elements = useMemo(
    () => (plot ? plot.getElementsByClassName('legendtoggle') : []),
    [plot]
  );

  const onHoverLegendToggle = useCallback(
    (e: any) => {
      // Re-find the element to grab the
      //   trace number
      for (var i = 0; i < elements.length; i++) {
        if (e.target !== elements[i]) continue;

        onHover(i);
        break;
      }
    },
    [onHover, elements]
  );

  const resetHoverLegendToggle = useCallback(() => {
    onHover(-1);
  }, [onHover]);

  const attachEventListeners = useCallback(() => {
    if (listenersAttached) return;

    for (var i = 0; i < elements.length; i++) {
      elements[i].addEventListener('mouseover', onHoverLegendToggle);
      elements[i].addEventListener('mouseout', resetHoverLegendToggle);
    }

    setListenersAttached(true);
  }, [
    elements,
    onHoverLegendToggle,
    resetHoverLegendToggle,
    listenersAttached
  ]);

  const removeEventListeners = useCallback(() => {
    for (var i = 0; i < elements.length; i++) {
      elements[i].removeEventListener('mouseover', onHoverLegendToggle);
      elements[i].removeEventListener('mouseout', resetHoverLegendToggle);
    }
  }, [elements, onHoverLegendToggle, resetHoverLegendToggle]);

  return {
    attachEventListeners,
    removeEventListeners
  };
}
