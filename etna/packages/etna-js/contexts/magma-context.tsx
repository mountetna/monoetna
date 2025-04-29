import React, { useContext, useState } from 'react';
import { ModelsObject } from '../models/magma-model';

export const MagmaContext = React.createContext({
  models: {} as ModelsObject,
  setModels: (models: ModelsObject) => {}
});

export const MagmaProvider = ({models: origModels, children}:{
  models?: ModelsObject;
  children: React.ReactNode;
}) => {
  const [ models, setModels ] = useState<ModelsObject>(origModels || {});

  return (<MagmaContext.Provider value={{
    models,
    setModels
  }}>
    {children}
  </MagmaContext.Provider>);
}
