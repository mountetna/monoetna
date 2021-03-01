import React, { useCallback, useState, useEffect } from 'react';
import ReactNotification, { store as notificationsStore } from 'react-notifications-component'
import 'react-notifications-component/dist/theme.css'
import 'animate.css/animate.min.css';

// Wrapper that provides the 'defaults' for etna ReactNotification container.  This wrapping does not do much now,
// but it does provide a common interface in which we can modify how all etna apps are configured.
export function Notifications() {
  return <ReactNotification />;
}

// Provides a set of functions that can add / remove notifications, accepting the options of
// react-notifications-component, while tracking 'local' notification state in such a way that
// a component using this hook
//   1.  Can freely use local 'identifiers' to add and remove notifications that do not conflict globally.
//   2.  Have all local notifications of the component using this hook automatically clear on dismount.
// This allows better 'scoping' of how notifications are shown.  One can still use this hook at a high level
// and 'share' the controller functions down to child components or in a React context if you wish to
// broaden the scope beyond just a single component's tree.
// https://github.com/teodosii/react-notifications-component#options
export function useNotifications() {
  const [localNotifications] = useState({});

  const removeAllLocalNotifications = useCallback(() => {
    Object.keys(localNotifications).forEach( k => {
      notificationsStore.removeNotification(localNotifications[k]);
      delete localNotifications[k];
    });
  }, []);

  const removeLocalNotification = useCallback((localId) => {
    const existing = localNotifications[localId];
    if (!existing) {
      return;
    }

    notificationsStore.removeNotification(existing);
    delete localNotifications[localId];
  }, [])

  const addLocalNotification = useCallback((localId, notificationOptions) => {
    const existing = localNotifications[localId];
    if (existing) {
      notificationsStore.removeNotification(existing);
    }

    localNotifications[localId] = notificationsStore.addNotification({
      ...notificationOptions
    });
  }, [])

  useEffect(() => {
    return removeAllLocalNotifications;
  }, []);

  return {
    removeAllLocalNotifications,
    addLocalNotification,
    removeLocalNotification,
  }
}
