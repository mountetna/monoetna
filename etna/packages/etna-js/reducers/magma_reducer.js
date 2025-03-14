/*
 *  models: {
 *
 *    model_1: {
 *      // this is a json document describing a Magma model
 *      template: {
 *      },
 *      documents: {
 *        // this is a json document describing a Magma record
 *        document_name1: {
 *          attribute1: value1
 *        }
 *      }
 *    },
 *    model_2: {
 *      etc.
 *    }
 *  }
 */
import {
  REVISE_DOCUMENT,
  DISCARD_REVISION,
  ADD_TEMPLATES_AND_DOCUMENTS,
  REMOVE_MODEL
} from '../actions/magma_actions';

import {UPLOAD_COMPLETE} from '../upload/workers/uploader';

var documents = function (old_documents, action) {
  if (!old_documents) old_documents = {};
  switch (action.type) {
    case ADD_TEMPLATES_AND_DOCUMENTS:
      let documents = {
        ...old_documents
      };
      for (let record_name in action.documents || {}) {
        documents[record_name] = {
          ...documents[record_name],
          ...action.documents[record_name]
        };
      }
      return documents;
    default:
      return old_documents;
  }
};

var revisions = function (old_revisions, action) {
  if (!old_revisions) old_revisions = {};
  switch (action.type) {
    case REVISE_DOCUMENT:
      return {
        ...old_revisions,
        [action.record_name]: {
          ...old_revisions[action.record_name],
          ...action.revision
        }
      };
    case DISCARD_REVISION:
      return {
        ...old_revisions,
        [action.record_name]: null
      };
    default:
      return old_revisions;
  }
};

var model = function (old_model, action) {
  if (!old_model)
    old_model = {
      documents: {},
      revisions: {},
      views: {}
    };

  switch (action.type) {
    case ADD_TEMPLATES_AND_DOCUMENTS:
      return {
        ...old_model,
        template: action.template ? action.template : old_model.template,
        documents: documents(old_model.documents, action)
      };
    case REVISE_DOCUMENT:
    case DISCARD_REVISION:
      return {
        ...old_model,
        revisions: revisions(old_model.revisions, action)
      };
    default:
      return old_model;
  }
};

var models = function (models, action) {
  if (!models) models = {};
  switch (action.type) {
    case ADD_TEMPLATES_AND_DOCUMENTS:
      return {
        ...models,
        ...Object.entries(action.models).reduce(
          (acc, [modelName, modelData]) => {
            acc[modelName] = model(models[modelName], {
              type: ADD_TEMPLATES_AND_DOCUMENTS,
              template: modelData.template,
              documents: modelData.documents
            });
            return acc;
          },
          {}
        )
      };
    case REMOVE_MODEL:
      let clone = {...models};
      delete clone[action.modelName];
      return clone;
    case REVISE_DOCUMENT:
    case DISCARD_REVISION:
      return {
        ...models,
        [action.model_name]: model(models[action.model_name], action)
      };
    default:
      return models;
  }
};

var magmaReducer = function (magma, action) {
  if (!magma)
    magma = {
      models: {},
      tables: {}
    };
  switch (action.type) {
    case ADD_TEMPLATES_AND_DOCUMENTS:
    case REVISE_DOCUMENT:
    case DISCARD_REVISION:
    case REMOVE_MODEL:
      return {
        ...magma,
        models: models(magma.models, action)
      };
    case UPLOAD_COMPLETE:
    // This case happens when a file upload completes.
    default:
      return magma;
  }
};

export default magmaReducer;
