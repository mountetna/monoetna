import YAML from 'yaml';

import Serializer from './interface';

export default class CwlSerializer extends Serializer {
  get json() {
    // Convert a CWL YAML file into JSON
    const doc = YAML.parseDocument(this.raw);
    return this.orderSteps(doc.toJSON());
  }

  orderSteps(json) {
    // Turn the json.steps into an Array with the steps
    //   in order, instead of a Hash. This is because the Hash
    //   keys don't come out of the YAML in guaranteed order.
    let orderedSteps = [];
    let allSteps = [];

    Object.keys(json.steps).forEach((stepName) => {
      allSteps.push(
        Object.assign({}, json.steps[stepName], {
          name: stepName
        })
      );
    });

    // We loop through allSteps and find all steps whose
    //   inputs are already in orderedSteps. Note that
    //   each step could have multiple input dependencies.
    // An input dependency is noted with a "/" in the
    //   input value. For example:
    //
    //   in: {
    //     data_from_previous_step: previous_step/output_variable_name
    //   }
    const numSteps = allSteps.length;
    while (orderedSteps.length < numSteps) {
      allSteps.forEach((step) => {
        if (orderedSteps.includes(step)) return;

        if (this.dependenciesSatisfied(step, orderedSteps))
          orderedSteps.push(step);
      });
    }

    json.steps = Object.assign([], orderedSteps);

    return json;
  }

  stepDependencies(step) {
    // Return a list of step names that the given
    //   step depends on for input. null otherwise.
    if ([] === step.in) return null;

    let dependencies = [];

    Object.keys(step.in).forEach((stepName) => {
      const inputVar = step.in[stepName];
      if (inputVar.includes('/')) dependencies.push(inputVar.split('/')[0]);
    });

    return dependencies.length > 0 ? dependencies : null;
  }

  dependenciesSatisfied(step, orderedSteps) {
    const dependencies = this.stepDependencies(step);

    if (null === dependencies) return true;

    const orderedStepNames = orderedSteps.map((s) => s.name);
    if (
      dependencies.every((dependency) => {
        return orderedStepNames.includes(dependency);
      })
    ) {
      return true;
    }

    return false;
  }
}
