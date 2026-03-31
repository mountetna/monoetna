import React, { useContext, useState, useCallback } from 'react';
import { ModelObject, ModelsObject } from '../models/magma-model';

export type CountsObject = {
  [modelName: string]: number;
}

export const MagmaContext = React.createContext({
  models: {} as ModelsObject,
  setModels: (models: ModelsObject) => {},
  setCounts: (counts: CountsObject) => {},
  counts: {} as CountsObject,
  mergeModels: (m1:ModelsObject, m2:ModelsObject):ModelsObject => ({} as ModelsObject)
});

export const MagmaProvider = ({models: origModels, children}:{
  models?: ModelsObject;
  children: React.ReactNode;
}) => {
  const [ models, setModels ] = useState<ModelsObject>(origModels || {});

  const mergeModel = (oldModel: ModelObject, newModel: ModelObject):ModelObject => ({
    ...oldModel,
    ...newModel,
    documents: {
      ...oldModel.documents,
      ...newModel.documents
    }
  });

  const mergeModels = (oldModels:ModelsObject,newModels:ModelsObject):ModelsObject => (
    {
      ...oldModels,
      ...Object.fromEntries(
        Object.keys(newModels).map(
          modelName => [ modelName, mergeModel(oldModels[modelName] || {}, newModels[modelName]) ]
        )
      )
    }
  );

  const [ counts, setCounts ] = useState<CountsObject>({} as CountsObject);

  return (<MagmaContext.Provider value={{
    models,
    setModels,
    mergeModels,
    counts,
    setCounts
  }}>
    {children}
  </MagmaContext.Provider>);
}
