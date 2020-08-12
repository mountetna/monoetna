import {useCallback} from 'react';
import {useNotifications} from "etna-js/components/Notifications";
import {useActionInvoker} from "etna-js/hooks/useActionInvoker";
import {requestDocuments} from "../actions/magma_actions";
import {documentOnWaiver, isWaiverPatientModel} from "../utils/patients";

export function useRequestDocuments() {
  const {
    removeAllLocalNotifications,
    addLocalNotification,
  } = useNotifications();

  const invoke = useActionInvoker();

  return useCallback((...args) => {
    removeAllLocalNotifications();

    return invoke(requestDocuments(...args)).then(payload => {
      const patient = Object.keys(payload.models).find(modelName => isWaiverPatientModel(modelName));
      if (patient) {
        const model = payload.models[patient];
        if (Object.values(model.documents).find(documentOnWaiver)) {
          addLocalNotification('waiver', {
            container: 'top-center',
            insert: 'bottom',
            type: 'warning',
            message: 'A patient returned in this dataset has not finished consent!',
            animationIn: ["animate__animated", "animate__fadeIn"],
            animationOut: ["animate__animated", "animate__fadeOut"],
            dismiss: false,
          })
        }
      }

      return payload;
    })
  }, [invoke]);
}