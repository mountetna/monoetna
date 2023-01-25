import React, {useState, createContext, useContext, useCallback} from 'react';
require('./ModalDialog.css');

const ModalDialogContext = createContext({});

export function ModalDialogContainer({children}) {
  const [modalState, setModalState] = useState({});
  const dismissModal = useCallback(() => setModalState({}), [setModalState]);
  const modalContextValue = {modalState, setModalState, dismissModal};
  const {component, opts} = modalState;

  function closeOnClickBackdrop() {
    return opts.hasOwnProperty('closeOnClickBackdrop')
      ? !!opts.closeOnClickBackdrop
      : true;
  }

  let modal = null;
  if (component) {
    console.log('opts', opts, closeOnClickBackdrop());
    modal = (
      <div
        className='etna-modal-window'
        onClick={() => (closeOnClickBackdrop() ? dismissModal() : null)}
      >
        <div onClick={(e) => e.stopPropagation()}>{component}</div>
      </div>
    );
  }

  return (
    <ModalDialogContext.Provider value={modalContextValue}>
      {children}
      {modal}
    </ModalDialogContext.Provider>
  );
}

export function useModal() {
  const {setModalState, dismissModal} = useContext(ModalDialogContext);
  const openModal = useCallback(
    function openModal(component, opts = {}) {
      setModalState({component, opts});
    },
    [setModalState]
  );

  return {openModal, dismissModal};
}
