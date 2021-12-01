import {showMessages} from './message_actions';
import {Exchange} from './exchange_actions';
import {
  getDocuments,
  getAnswer,
  getTSVForm,
  getQueryTSVForm,
  postRevisions
} from '../api/magma_api';
import {setMetisCookie} from './metis_actions';
import {filePathComponents} from '../selectors/magma';
import {TEMP} from './file_actions';
import {dispatchUploadWork} from '../upload/actions/upload_actions';
import {AddUploadCommand, Upload} from '../upload/workers/uploader';

export const REVISE_DOCUMENT = 'REVISE_DOCUMENT';
export const DISCARD_REVISION = 'DISCARD_REVISION';
export const ADD_PREDICATES = 'ADD_PREDICATES;';
export const ADD_TEMPLATES_AND_DOCUMENTS = 'ADD_TEMPLATES_AND_DOCUMENTS';

export const reviseDocument = (
  document,
  template,
  attribute,
  revised_value
) => {
  return {
    type: REVISE_DOCUMENT,
    model_name: template.name,
    record_name: document[template.identifier],
    revision: {
      [attribute.name]: revised_value
    }
  };
};

export const discardRevision = (record_name, model_name) => {
  return {
    type: DISCARD_REVISION,
    model_name: model_name,
    record_name: record_name
  };
};

export const addPredicates = (predicates) => {
  return {
    type: ADD_PREDICATES,
    predicates
  };
};

export const addTemplatesAndDocuments = (models) => {
  return {
    type: ADD_TEMPLATES_AND_DOCUMENTS,
    models
  }
}

/*
 * Here we add the models and documents to the store. At the same time we strip
 * off the model namespacing. The server returns the full name of the model.
 * And we want that behavior so we can know from which project the data is
 * coming from. However, we do not want to propigate that namespacing to the UI.
 */
export const consumePayload = (dispatch, response) => {
  if (response.models) {
    dispatch(addTemplatesAndDocuments(response.models));
  }
};

const showError = dispatch => response => {
  if ('error' in response) {
    dispatch(showMessages([`There was a ${response.type} error.`]));
    console.log(response.error);
  }
}

export const requestDocuments = ({
  model_name,
  record_names,
  attribute_names,
  filter,
  show_disconnected,
  page,
  page_size,
  collapse_tables,
  exchange_name,
  output_predicate
}) => dispatch => {
  const exchange = new Exchange(dispatch, exchange_name);
  return getDocuments(
    {
      model_name,
      record_names,
      attribute_names,
      filter,
      show_disconnected,
      page,
      page_size,
      collapse_tables,
      output_predicate
    },
    exchange.fetch.bind(exchange)
  ).then( response => {
    showError(dispatch)(response);
    consumePayload(dispatch, response);
    return Promise.resolve(response);
  }).catch(e => {
    return Promise.resolve(e).then((response) => {
      let errStr = response.error
        ? response.error
        : response.errors
        ? response.errors.map((error) => `* ${error}`)
        : response;
      errStr = [`### Our request was refused.\n\n${errStr}`];
      dispatch(showMessages(errStr));
      return Promise.reject(e);
    });
  })
}

export const requestModels = () => requestDocuments({
  model_name: 'all',
  record_names: [],
  attribute_names: 'all',
  exchange_name: 'request-models'
});

export const requestModel = (model_name) => requestDocuments({
  model_name: model_name,
  record_names: [],
  attribute_names: 'all',
  exchange_name: 'request-model'
});

export const requestIdentifiers = () => requestDocuments({
  model_name: 'all',
  record_names: 'all',
  attribute_names: 'identifier',
  exchange_name: 'request-identifiers'
});

const uploadFileRevisions = (model_name, revisions, response, dispatch) => {
  // Handle any file upload revisions by sending the file to Metis.
  // Here Magma's response should be a temporary upload URL.
  // 1. Find the file attributes from the response.models.model_name.template.attributes.
  // 2. Grab the file(s) from the original revisions list for each record.
  // 3. Upload the file(s) to Metis.
  // 4. Send another update to Magma with the file location(s) and
  //      original filenames.
  // 5. Update the response with final Magma URLs for each file attribute.
  // 6. Return the response via a Promise.
  uploadTemporaryFiles(model_name, revisions, response, dispatch);

  return new Promise((resolve, reject) => {
    resolve(response);
  });
};

const attributeNamesByTypeForModel = (modelAttributes, attribute_type) => {
  const fileAttributeNames = [];
  Object.keys(modelAttributes).forEach((attribute_name) => {
    let attribute = modelAttributes[attribute_name];
    if (attribute.attribute_type === attribute_type) {
      fileAttributeNames.push(attribute.attribute_name);
    }
  });
  return fileAttributeNames;
};

const getCookie = (cookieName) => {
  const cookies = decodeURIComponent(document.cookie).split(';');
  for (var i = 0; i < cookies.length; i++) {
    let cookie = cookies[i].split('=');
    if (cookie[0].trim() === cookieName) {
      return cookie[1];
    }
  }
  return null;
};

const uploadTemporaryFiles = (model_name, revisions, response, dispatch) => {
  const model = response.models[model_name];
  if (!model) return;

  const modelAttributes = model.template.attributes;
  const fileAttributeNames = attributeNamesByTypeForModel(
    modelAttributes,
    'file'
  ).concat(attributeNamesByTypeForModel(modelAttributes, 'image'));

  Object.keys(revisions).forEach((record_name) => {
    const revision = revisions[record_name];
    Object.keys(revision).forEach((attribute_name) => {
      if (fileAttributeNames.indexOf(attribute_name) === -1) {
        return;
      }

      const fileRevision = revision[attribute_name];

      if (fileRevision.path !== TEMP) {
        return;
      }

      const tempRevision =
        response.models[model_name].documents[record_name][attribute_name];

      let {hostname, file_name} = filePathComponents(tempRevision.path);

      // To avoid collisions, we need to create a new file with name
      //   matching the random name that Magma returns.
      const originalFile = fileRevision.original_files[0];
      const tempFile = new File([originalFile], file_name, {
        type: originalFile.type
      });

      setMetisCookie(dispatch, `https://${hostname}`).then(() => {
        dispatchUploadWork(
          dispatch,
          AddUploadCommand(
            Upload({
              file_name,
              url: tempRevision.path,
              file: tempFile,
              project_name: CONFIG.project_name,
              metis_uid: getCookie(CONFIG.metis_uid_name),
              attribute_name,
              model_name,
              record_name,
              original_filename: originalFile.name
            })
          )
        );
      });
    });
  });

  return;
};

export const finalizeUpload = (
  model_name,
  template,
  record_name,
  attribute_name,
  upload
) => {
  return (dispatch) => {
    let {project_name, bucket_name, file_name} = filePathComponents(upload.url);

    sendRevisions(model_name, template, {
      [record_name]: {
        [attribute_name]: {
          path: `metis://${project_name}/${bucket_name}/${file_name}`,
          original_filename: upload.original_filename
        }
      }
    })(dispatch);
  };
};

export const sendRevisions = (
  model_name,
  model_template,
  revisions,
  success,
  error
) => dispatch => {
  const exchange =  new Exchange(dispatch, `revisions-${model_name}`);
  return postRevisions(
    formatRevisions(revisions, model_name, model_template),
    exchange.fetch.bind(exchange)
  ).then(
    response => uploadFileRevisions(model_name, revisions, response, dispatch).then(
      (updatedResponse) => {
        consumePayload(dispatch, updatedResponse);
        for (var record_name in revisions) {
          dispatch(discardRevision(record_name, model_name));
        }

        if (success != undefined) success();
      }
    )
  ).catch(e => {
    e.then((response) => {
      let errStr = response.error
        ? response.error
        : response.errors.map((error) => `* ${error}`);
      errStr = [`### The change we sought did not occur.\n\n${errStr}`];
      dispatch(showMessages(errStr));
    });

    if (error != undefined) error();
  });
}

// export for testing
export const formatRevisions = (revisions, model_name, model_template) => {
  const modelRevisions = {};
  modelRevisions[model_name] = cleanFileCollectionRevisions(
    revisions,
    model_template
  );

  const formattedRevs = {
    revisions: modelRevisions
  };
  return formattedRevs;
};

const cleanFileCollectionRevisions = (revisions, model_template) => {
  // Need to strip out any empty string FileCollection attribute updates.
  // If someone clicks + but does not select a file or set a Metis path,
  //    there will be a "" value in the FileCollection revision array.
  // Magma will throw a ServerError on that value, because strings cannot
  //    be set for FileCollection values.
  const modelAttributes = model_template.attributes;
  const fileCollectionAttributeNames = attributeNamesByTypeForModel(
    modelAttributes,
    'file_collection'
  );

  let cleanRevisions = {};

  Object.keys(revisions).forEach((record_name) => {
    const revision = revisions[record_name];
    cleanRevisions[record_name] = {};
    Object.keys(revision).forEach((attribute_name) => {
      if (fileCollectionAttributeNames.indexOf(attribute_name) === -1) {
        cleanRevisions[record_name][attribute_name] = revision[attribute_name];
        return;
      }

      // Remove empty strings
      cleanRevisions[record_name][attribute_name] = revision[
        attribute_name
      ].filter((rev) => {
        return '' !== rev;
      });
    });
  });

  return cleanRevisions;
};

// Download a TSV from magma via Timur via /retrieve
export const requestTSV = (params) => (dispatch) => getTSVForm(params);

export const requestQueryTSV = (params) => (dispatch) => getQueryTSVForm(params);

export const requestAnswer = (question, callback = null) => dispatch => {
  const exchange = new Exchange(
    dispatch,
    Array.isArray(question) ? JSON.stringify(question) : question
  );

  return getAnswer(question, exchange.fetch.bind(exchange)).then(
    response => {
      showError(dispatch)(response);
      if (callback) callback(response);
      else return response;
  }).catch( error => {
    console.log(error);
    throw new Error(error);
  });
}

export const requestPredicates = () => {
  return (dispatch) => {
    let localCallback = (response) => {
      dispatch(addPredicates(response.predicates));
    };

    dispatch(requestAnswer('::predicates', localCallback));
  };
};
