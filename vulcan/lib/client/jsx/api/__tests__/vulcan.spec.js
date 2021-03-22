const fs = require('fs');
const path = require('path');

import nock from 'nock';

import {getWorkflows, submit, getData, downloadUrlUpdated} from '../vulcan';
import {mockStore, stubUrl, mockFetch, cleanStubs} from 'etna-js/spec/helpers';

describe('Vulcan API', () => {
  afterEach(() => {
    cleanStubs();
    nock.cleanAll();
  });
  afterAll(nock.restore);

  mockFetch();

  it('getWorkflows returns all workflows', (done) => {
    const expectedWorkflows = {
      analysis_1: {},
      analysis_2: {},
      analysis_3: {}
    };

    stubUrl({
      verb: 'get',
      path: ROUTES.fetch_workflows(),
      response: expectedWorkflows,
      headers: {
        'Content-type': 'application/json'
      },
      host: CONFIG.vulcan_host
    });

    return getWorkflows().then((data) => {
      expect(data).toEqual(expectedWorkflows);
      done();
    });
  });

  xit('getWorkflow returns a YAML workflow by default', (done) => {
    const sample_yaml = fs.readFileSync(
      path.resolve(__dirname, '../../spec/fixtures/sample_cwl.yaml'),
      'utf8'
    );

    stubUrl({
      verb: 'get',
      path: ROUTES.fetch_workflow('test'),
      response: sample_yaml,
      headers: {
        'Content-type': 'text/yaml'
      },
      host: CONFIG.vulcan_host
    });

    return getWorkflow('test').then((data) => {
      expect(data).toEqual(sample_yaml);
      done();
    });
  });

  xit('getWorkflow returns a JSON workflow when requested', (done) => {
    const workflow = {
      class: 'Workflow',
      cwlVersion: 'v1.1',
      inputs: {
        bool_input: {
          default: true,
          label: 'Sample boolean',
          type: 'boolean'
        },
        int_input: {
          default: 42,
          label: 'Sample input',
          type: 'int'
        }
      },
      outputs: {
        sample_data: {
          outputSource: 'final_step/sample_data',
          type: 'File'
        }
      },
      steps: [
        [
          {
            in: [],
            out: ['choice_set'],
            run: 'first_step.cwl',
            name: 'first_step'
          },
          {
            in: {
              all_pool_names: 'first_step/choice_set'
            },
            out: ['subset'],
            run: 'ui_pick_subset.cwl',
            name: 'ui_pick_subset'
          },
          {
            in: {
              bool_input: 'bool_input',
              data: 'ui_pick_subset/subset',
              int_input: 'int_input'
            },
            out: ['sample_data'],
            run: 'final_step.cwl',
            name: 'final_step'
          }
        ]
      ]
    };

    const json = true;

    stubUrl({
      verb: 'get',
      path: ROUTES.fetch_workflow('test', json),
      response: workflow,
      headers: {
        'Content-type': 'application/json'
      },
      host: CONFIG.vulcan_host
    });

    return getWorkflow('test', json).then((data) => {
      expect(data).toEqual(workflow);
      done();
    });
  });

  it('submit posts existing inputs and fetches getData', (done) => {
    const inputs = [
      [
        {
          name: 'first_step'
        },
        {
          all_pool_names: [1, 2, 3, 4],
          name: 'ui_pick_subset'
        },
        {
          name: 'final_step'
        }
      ],
      [
        {
          name: 'first_step'
        },
        {
          all_pool_names: [1],
          name: 'branching_step'
        },
        {
          name: 'final_step'
        }
      ]
    ];

    const url = '/api/workflows/test/file1.txt';
    const data = [1, 2, 4, 'abc'];
    const status = [
      [
        {
          name: 'first_step',
          status: 'complete',
          downloads: {
            sum: `${CONFIG.vulcan_host}${url}`
          }
        },
        {
          all_pool_names: [1, 2, 3, 4],
          name: 'ui_pick_subset',
          status: 'pending'
        },
        {
          name: 'final_step',
          status: 'pending'
        }
      ],
      [
        {
          name: 'first_step',
          status: 'complete',
          downloads: {
            sum: `${CONFIG.vulcan_host}${url}`
          }
        },
        {
          all_pool_names: [1],
          name: 'branching_step',
          status: 'complete'
        },
        {
          name: 'final_step',
          status: 'pending'
        }
      ]
    ];

    stubUrl({
      verb: 'post',
      path: ROUTES.submit('example', 'test'),
      request: () => ({inputs, key: 'session_key'}),
      response: {
        status,
        session: {}
      },
      headers: {
        'Content-type': 'application/json'
      },
      host: CONFIG.vulcan_host
    });
    stubUrl({
      verb: 'get',
      path: url,
      response: data,
      headers: {
        'Content-type': 'application/json'
      },
      host: CONFIG.vulcan_host
    });

    const mockSetSession = jest.fn();
    const mockSetStatus = jest.fn();
    const mockSetData = jest.fn();

    let context = {
      workflow: {
        name: 'test',
        steps: [
          [
            {
              name: 'first_step',
              out: ['data']
            },
            {
              name: 'ui_pick_subset',
              run: 'ui-queries/select.cwl',
              in: [{source: ['first_step', 'data']}]
            },
            {
              name: 'final_step'
            }
          ]
        ]
      },
      session: {
        inputs,
        key: 'session_key'
      },
      status: [
        [{name: 'first_step'}, {name: 'ui_pick_subset'}, {name: 'final_step'}]
      ],
      pathIndex: 0,
      setSession: mockSetSession,
      setStatus: mockSetStatus,
      setData: mockSetData
    };

    return submit(context).then(() => {
      expect(mockSetSession.mock.calls.length).toBe(1);
      expect(mockSetStatus.mock.calls.length).toBe(1);
      expect(mockSetData.mock.calls.length).toBe(1);
      done();
    });
  });

  it('submit does not refetch data if output URL stays same', (done) => {
    const inputs = [
      [
        {
          name: 'first_step'
        }
      ]
    ];

    const url = '/api/workflows/test/file1.txt';
    const data = [1, 2, 4, 'abc'];
    const status1 = [
      [
        {
          name: 'first_step',
          status: 'complete',
          downloads: {
            sum: `${CONFIG.vulcan_host}${url}`
          },
          data: {
            sum: data
          }
        }
      ]
    ];

    stubUrl({
      verb: 'post',
      path: ROUTES.submit('example', 'test'),
      request: () => ({inputs, key: 'session_key'}),
      response: {
        status: status1,
        session: {}
      },
      headers: {
        'Content-type': 'application/json'
      },
      host: CONFIG.vulcan_host
    });

    const mockSetSession = jest.fn();
    const mockSetStatus = jest.fn();
    const mockSetData = jest.fn();

    let context = {
      workflow: {
        name: 'test',
        steps: [
          [
            {name: 'first_step', out: ['data']},
            {
              name: 'second_step',
              run: 'ui-queries/select.cwl',
              in: [{source: ['first_step', 'data']}]
            }
          ]
        ]
      },
      session: {
        inputs,
        key: 'session_key'
      },
      status: status1,
      pathIndex: 0,
      setSession: mockSetSession,
      setStatus: mockSetStatus,
      setData: mockSetData
    };

    return submit(context).then(() => {
      expect(mockSetSession.mock.calls.length).toBe(1);
      expect(mockSetStatus.mock.calls.length).toBe(1);
      expect(mockSetData.mock.calls.length).toBe(0);
      done();
    });
  });

  it('submit refetches data if output URL changes', (done) => {
    const inputs = [
      [
        {
          name: 'first_step'
        }
      ]
    ];

    const url = '/api/workflows/test/file1.txt';
    const url2 = '/api/workflows/test/file2.txt';
    const data = [1, 2, 4, 'abc'];
    const data2 = [5, 6, 7, 'xyz'];
    const status1 = [
      [
        {
          name: 'first_step',
          status: 'complete',
          downloads: {
            sum: `${CONFIG.vulcan_host}${url}`
          },
          data: {
            sum: data
          }
        }
      ]
    ];
    const status2 = [
      [
        {
          name: 'first_step',
          status: 'complete',
          downloads: {
            sum: `${CONFIG.vulcan_host}${url2}`
          }
        }
      ]
    ];

    stubUrl({
      verb: 'post',
      path: ROUTES.submit('example', 'test'),
      request: () => ({inputs, key: 'session_key'}),
      response: {
        status: status2,
        session: {}
      },
      headers: {
        'Content-type': 'application/json'
      },
      host: CONFIG.vulcan_host
    });
    stubUrl({
      verb: 'get',
      path: url2,
      response: data2,
      headers: {
        'Content-type': 'application/json'
      },
      host: CONFIG.vulcan_host
    });

    const mockSetSession = jest.fn();
    const mockSetStatus = jest.fn();
    const mockSetData = jest.fn();

    let context = {
      workflow: {
        name: 'test',
        steps: [
          [
            {name: 'first_step', out: ['data']},
            {
              name: 'second_step',
              run: 'ui-queries/select.cwl',
              in: [{source: ['first_step', 'data']}]
            }
          ]
        ]
      },
      session: {
        inputs,
        key: 'session_key'
      },
      status: status1,
      pathIndex: 0,
      setSession: mockSetSession,
      setStatus: mockSetStatus,
      setData: mockSetData
    };

    return submit(context).then(() => {
      expect(mockSetSession.mock.calls.length).toBe(1);
      expect(mockSetStatus.mock.calls.length).toBe(1);
      expect(mockSetData.mock.calls.length).toBe(1);
      done();
    });
  });

  it('submit does not fetch data for non-UI related steps', (done) => {
    const inputs = [
      [
        {
          name: 'first_step'
        }
      ]
    ];

    const url = '/api/workflows/test/file1.txt';
    const url2 = '/api/workflows/test/file2.txt';
    const data = [1, 2, 4, 'abc'];
    const data2 = [5, 6, 7, 'xyz'];
    const status1 = [
      [
        {
          name: 'first_step',
          status: 'complete',
          downloads: {
            sum: `${CONFIG.vulcan_host}${url}`
          },
          data: {
            sum: data
          }
        }
      ]
    ];
    const status2 = [
      [
        {
          name: 'first_step',
          status: 'complete',
          downloads: {
            sum: `${CONFIG.vulcan_host}${url2}`
          }
        }
      ]
    ];

    stubUrl({
      verb: 'post',
      path: ROUTES.submit('example', 'test'),
      request: () => ({inputs, key: 'session_key'}),
      response: {
        status: status2,
        session: {}
      },
      headers: {
        'Content-type': 'application/json'
      },
      host: CONFIG.vulcan_host
    });
    stubUrl({
      verb: 'get',
      path: url2,
      response: data2,
      headers: {
        'Content-type': 'application/json'
      },
      host: CONFIG.vulcan_host
    });

    const mockSetSession = jest.fn();
    const mockSetStatus = jest.fn();
    const mockSetData = jest.fn();

    let context = {
      workflow: {
        name: 'test',
        steps: [
          [
            {name: 'first_step', out: ['data']},
            {
              name: 'second_step',
              run: 'scripts/select.cwl',
              in: [{source: ['first_step', 'data']}]
            }
          ]
        ]
      },
      session: {
        inputs,
        key: 'session_key'
      },
      status: status1,
      pathIndex: 0,
      setSession: mockSetSession,
      setStatus: mockSetStatus,
      setData: mockSetData
    };

    return submit(context).then(() => {
      expect(mockSetSession.mock.calls.length).toBe(1);
      expect(mockSetStatus.mock.calls.length).toBe(1);
      expect(mockSetData.mock.calls.length).toBe(0);
      done();
    });
  });

  it('getData returns data from the specified location', (done) => {
    const url = '/api/workflows/test/file1.txt';
    const data = [1, 2, 4, 'abc'];

    stubUrl({
      verb: 'get',
      path: url,
      response: data,
      headers: {
        'Content-type': 'application/json'
      },
      host: CONFIG.vulcan_host
    });

    return getData(`${CONFIG.vulcan_host}${url}`).then((returnedData) => {
      expect(data).toEqual(returnedData);
      done();
    });
  });

  describe('downloadUrlUpdated', () => {
    it('returns true if data does not exist on new Status', () => {
      let oldStatus = {
        name: 'foo'
      };

      let newStatus = {
        downloads: {
          key: 'URL'
        }
      };

      let result = downloadUrlUpdated(oldStatus, newStatus, 'key');
      expect(result).toEqual(true);

      oldStatus = {
        name: 'foo',
        data: {
          another_key: 'blob'
        }
      };

      result = downloadUrlUpdated(oldStatus, newStatus, 'key');
      expect(result).toEqual(true);
    });

    it('returns false if url did not change', () => {
      let oldStatus = {
        name: 'foo',
        downloads: {
          key: 'URL'
        },
        data: {
          key: 'blob'
        }
      };

      let newStatus = {
        downloads: {
          key: 'URL'
        },
        data: {
          key: 'blob'
        }
      };

      let result = downloadUrlUpdated(oldStatus, newStatus, 'key');
      expect(result).toEqual(false);
    });

    it('returns true if url changed but data exists', () => {
      let oldStatus = {
        name: 'foo',
        downloads: {
          key: 'URL'
        },
        data: {
          key: 'blob'
        }
      };

      let newStatus = {
        downloads: {
          key: 'URL2'
        },
        data: {
          key: 'blob2'
        }
      };

      let result = downloadUrlUpdated(oldStatus, newStatus, 'key');
      expect(result).toEqual(true);
    });

    it('returns true if no downloads for the old step', () => {
      let oldStatus = {
        name: 'foo'
      };

      let newStatus = {
        downloads: {
          key: 'URL2'
        },
        data: {
          key: 'blob2'
        }
      };

      let result = downloadUrlUpdated(oldStatus, newStatus, 'key');
      expect(result).toEqual(true);
    });
  });
});
