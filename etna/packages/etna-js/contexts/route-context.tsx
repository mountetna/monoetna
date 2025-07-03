import React, { useContext, useState } from 'react';

export const RouteContext = React.createContext({});

export const RouteProvider = ({children, state}:{
  state: any;
  children: React.ReactNode;
}) => {
  return (<RouteContext.Provider value={{
  }}>
    {children}
  </RouteContext.Provider>);
}
