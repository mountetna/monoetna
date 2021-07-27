import {Cancellable} from "etna-js/utils/cancellable";
import {useCallback} from "react";

export const defaultConfirmationHelpers = {
  confirm(message: string, context: Cancellable) {
    return Promise.resolve(true)
  },
}

export function useConfirmation(windowConfirm = window.confirm): typeof defaultConfirmationHelpers {
  const openConfirmationModal = useCallback((message: string): [Promise<boolean>, Function] => {
    // TODO: Implement some proper modal that allows for closing via the returned callback.
    return [Promise.resolve(windowConfirm(message)), () => {
      // This callback should close the modal if displaying a truly async modal.
      // Used to cleanup when the modal needs to be forced closed.
    }]
  }, [windowConfirm]);

  const confirm = useCallback((message: string, context: Cancellable) => {
    const [promise, cleanup] = openConfirmationModal(message);
    return context.race(promise).then(({ cancelled, result }) => {
      if (cancelled || typeof result === "undefined") {
        cleanup();
        return false;
      }

      return result;
    })
  }, [openConfirmationModal])

  return {
    confirm,
  }
}