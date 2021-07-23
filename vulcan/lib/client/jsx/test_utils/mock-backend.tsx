import {
  defaultSessionStatusResponse, defaultStepStatus, SessionStatusResponse, VulcanSession, Workflow, WorkflowStep,
} from "../api_types";
import {ProviderProps, VulcanContextData} from "../contexts/vulcan_context";
import {
  sourceNameOfReference, splitSource, stepOfStatus, uiOutputOfStep, uiQueryOfStep
} from "../selectors/workflow_selectors";
import {DataEnvelope} from "../components/workflow/user_interactions/inputs/input_types";
import objectHash from "object-hash";
import {mapSome, Maybe, some} from "../selectors/maybe";

export type Overrides = Partial<VulcanContextData> & Partial<ProviderProps>;
export type Script = (inputs: DataEnvelope<any>) => Promise<DataEnvelope<any | ArrayBuffer>>;

function asJsonOrBlob(ab: ArrayBuffer): any | Blob {
  try {
    return JSON.parse(new TextDecoder().decode(ab));
  } catch (e) {
    return new Blob([ab]);
  }
}

function asArrayBuffer(json: any): ArrayBuffer {
  if (json instanceof ArrayBuffer) return json;
  const string = JSON.stringify(json)
  return new TextEncoder().encode(string).buffer;
}

export function createFakeBackend(workflows: Workflow[], scripts: DataEnvelope<Script>,): Overrides {
  const dataByUrl: DataEnvelope<ArrayBuffer> = {};
  const completed: DataEnvelope<MockBuildTarget> = {};
  const errored: DataEnvelope<string> = {};
  const pending: DataEnvelope<boolean> = {};

  return {
    async getWorkflows() {
      return {workflows};
    }, async getData(url: string): Promise<any> {
      if (url in dataByUrl) {
        return asJsonOrBlob(dataByUrl[url]);
      }

      throw new Error(`Could not find data for url ${url}`)
    }, async pollStatus(session: VulcanSession): Promise<SessionStatusResponse> {
      const workflow = getWorkflow(session.workflow_name);
      const response = {...defaultSessionStatusResponse, session};
      return updateStatus(workflow, session, response)
    }, async postInputs(session: VulcanSession): Promise<SessionStatusResponse> {
      const workflow = getWorkflow(session.workflow_name);
      const response = {...defaultSessionStatusResponse, session};
      scheduleNewWork(workflow, session);
      return updateStatus(workflow, session, response);
    }
  }

  function gatherInputs(step: WorkflowStep, session: VulcanSession, workflow: Workflow): Maybe<DataEnvelope<any>> {
    const inputs: DataEnvelope<any> = {};
    for (let input of step.in) {
      const {source, id} = input;

      if (source in session.inputs) {
        inputs[id] = session.inputs[source];
        continue;
      }

      const [sourceStepName, sourceOutputName] = splitSource(source);
      const bt = sourceStepName ?
        getBuildTargetFor(sourceStepName, session, workflow) :
        getPrimaryInputsBuildTarget(workflow, session);

      const hash = bt.hash();
      if (hash in completed) {
        inputs[id] = asJsonOrBlob(dataByUrl[urlOf(hash, sourceOutputName)])
      } else {
        // Cannot be run, not all inputs satisfied.
        return null;
      }
    }

    return some(inputs);
  }

  function scheduleNewWork(workflow: Workflow, session: VulcanSession) {
    // Try running any step that can be run.
    workflow.steps[0].forEach(step => {
      console.log('trying to scheduling', step.name);
      if (uiQueryOfStep(step)) return;
      if (uiOutputOfStep(step)) return;

      const bt = getBuildTargetFor(step.name, session, workflow);
      const hash = bt.hash();

      if (!(hash in pending) && !(hash in completed)) {
        const inputs = gatherInputs(step, session, workflow);

        mapSome(inputs, inputs => {
          // Clear the last error
          errored[hash] = "";
          pending[hash] = true;

          bt.run(inputs).finally(() => {
            pending[hash] = false;
          }).then(results => {
            Object.entries(results).forEach(([outputName, data]) => {
              dataByUrl[urlOf(hash, outputName)] = data;
            })
            completed[hash] = bt;
            scheduleNewWork(workflow, session);
          }).catch(e => {
            errored[hash] = e + ""
          })
        })
      }
    })
  }


  function updateStatus(workflow: Workflow, session: VulcanSession, response: SessionStatusResponse) {
    response.status = [
      workflow.steps[0].map(step => {
        const bt = getBuildTargetFor(step.name, session, workflow);
        const hash = bt.hash();

        const status = {...defaultStepStatus};
        status.hash = hash;
        status.name = step.name;

        if (uiOutputOfStep(step) || uiQueryOfStep(step)) {
          if (gatherInputs(step, session, workflow) != null) {
            status.status = 'complete';
            return status;
          }
        }

        const error = errored[hash];
        const complete = completed[hash];
        const running = pending[hash];

        if (running) {
          status.status = 'running';
        } else if (error) {
          status.status = 'error';
          status.error = error;
        } else if (complete) {
          status.status = 'complete';
          const downloads = status.downloads = {} as DataEnvelope<string>;
          step.out.forEach(outputName => (downloads[outputName] = urlOf(hash, outputName)));
        } else {
          status.status = 'pending';
        }

        return status;
      })
    ];

    return response;
  }

  function getPrimaryInputsBuildTarget(workflow: Workflow, session: VulcanSession) {
    const inputs: MockStorageFile[] = [];
    Object.keys(workflow.inputs).forEach(primaryInput => {
      if (primaryInput in session.inputs) {
        inputs.push(new MockStorageFile(
          primaryInput,
          primaryInput,
          objectHash.MD5({fulfilled: session.inputs[primaryInput]})
        ))
      } else {
        inputs.push(new MockStorageFile(primaryInput, primaryInput, objectHash.MD5(null)))
      }
    })

    return new MockBuildTarget(Object.keys(workflow.inputs), inputs, async () => ({}));
  }

  function getBuildTargetFor(stepName: string,
    session: VulcanSession,
    workflow: Workflow,
    cache: DataEnvelope<MockBuildTarget> = {},
  ): MockBuildTarget {
    const step = stepOfStatus(stepName, workflow);
    if (!step) throw new Error(`Step ${stepName} could not be resolved in workflow ${workflow.name}`)
    if (stepName in cache) return cache[stepName];

    const inputs: MockStorageFile[] = [];

    const uiStep = !!uiOutputOfStep(step) || !!(uiQueryOfStep(step));

    if (uiStep) {
      for (let output of step.out) {
        const source = sourceNameOfReference([step.name, output]);
        if (source in session.inputs) {
          inputs.push(new MockStorageFile(output, output, objectHash.MD5({fulfilled: session.inputs[source]})));
        } else {
          inputs.push(new MockStorageFile(output, output, objectHash.MD5(null)));
        }
      }
    }

    for (let input of step.in) {
      const {source, id} = input;
      const [stepName, outputName] = splitSource(source);
      const bt = stepName ?
        getBuildTargetFor(stepName, session, workflow, cache) :
        getPrimaryInputsBuildTarget(workflow, session);

      inputs.push(bt.outputs[outputName].asInput(id));
    }

    return (cache[step.name] = new MockBuildTarget(step.out, inputs, uiStep ? async () => ({}) : getScript(step.run),))
  }

  function getWorkflow(name: string) {
    const workflow = workflows.find(({name: n}) => name === n);
    if (!workflow) throw new Error(`Workflow ${name} does not exist in the mock!`);
    return workflow;
  }

  function getScript(run: string) {
    const script = scripts[run];
    if (!script) throw new Error(`${run} does not exist in the mock!`);
    return script;
  }

  function urlOf(hash: string, outputName: string) {
    return `http://data/${hash}/${outputName}`
  }
}


interface Hashable {
  hash(): string;
}

class MockBuildTarget {
  constructor(private outputFileNames: string[], private inputFiles: MockStorageFile[], private script: Script,) {
  }

  hash() {
    const {outputFileNames, inputFiles, script} = this;
    const scriptAsString = script.toString();
    const inputFilesAsStrings = inputFiles.map(f => f.hash());

    return objectHash.MD5({
      outputFileNames, scriptAsString, inputFilesAsStrings,
    })
  }

  get outputs() {
    const outputs: DataEnvelope<MockStorageFile> = {};
    this.outputFileNames.forEach(fileName => {
      outputs[fileName] = new MockStorageFile(fileName, fileName, this.hash());
    })

    return outputs;
  }

  async run(inputs: DataEnvelope<any>) {
    const outputs: DataEnvelope<ArrayBuffer> = {};
    const result = await this.script(inputs);
    this.outputFileNames.forEach(fileName => {
      if (fileName in result) {
        let value = result[fileName];

        let ab: ArrayBuffer = value;
        if (!(value instanceof ArrayBuffer)) {
          const str = JSON.stringify(value);
          ab = new TextEncoder().encode(str).buffer;
        }

        outputs[fileName] = ab;
      }
    })

    return outputs;
  }
}

class MockStorageFile implements Hashable {
  constructor(public sourceName: string, public inputName: string, public targetHash: string) {
  }

  asInput(inputName: string) {
    return new MockStorageFile(this.sourceName, inputName, this.targetHash);
  }

  hash() {
    const {sourceName, inputName, targetHash} = this;
    return objectHash.MD5({sourceName, inputName, targetHash});
  }
}

//
// function m() {
//   const bytes = new Uint8Array(contents);
//   let binary: string[] = [];
//   const len = bytes.byteLength;
//   for (let i = 0; i < len; i++) {
//     binary.push(String.fromCharCode(bytes[i]));
//   }
//   const contentsString = btoa(binary.join(''));
//   return objectHash.MD5({sourceName, inputName, contentsString});
// }