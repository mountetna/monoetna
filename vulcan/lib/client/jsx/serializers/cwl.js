import YAML from 'yaml';

import Serializer from './interface';

export default class CwlSerializer extends Serializer {
  get json() {
    // Convert a CWL YAML file into JSON
    const doc = YAML.parseDocument(this.raw);
    return injectStepMetadata(doc.toJSON());
  }

  injectStepMetadata(json) {
    // Add some logic to add stepNumber to each Step.
    // These have to be calculated by graph traversal to figure
    //   out which inputs are dependent on which other steps.
    let orderedSteps = [];
    let allSteps = [];

    Object.keys(json.steps).forEach((stepName) => {
      allSteps.push(
        Object.assign({}, json.steps[stepName], {
          name: stepName
        })
      );
    });

    while (orderedSteps.length < allSteps.length) {}

    orderedSteps.forEach((step) => {
      json.steps[step.name] = Object.assign({}, step);
    });

    return json;
  }
}
