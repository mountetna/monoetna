import {
  completedUiStepsSelector,
  nextUiStepsSelector,
  errorStepsSelector
} from '../workflow';

describe('Workflow Selectors', () => {
  describe('completedUiStepsSelector', () => {
    it('returns only completed UI steps', () => {
      const context = {
        status: [
          [
            {
              name: 'step1',
              status: 'pending'
            },
            {
              name: 'step2',
              status: 'complete'
            },
            {
              name: 'step3',
              status: 'complete'
            },
            {
              name: 'step4',
              status: 'error'
            }
          ]
        ],
        pathIndex: 0,
        workflow: {
          steps: [
            [
              {name: 'step1', run: 'ui-queries/something.cwl'},
              {name: 'step2', run: 'scripts/real.cwl'},
              {name: 'step3', run: 'ui-queries/other.cwl'},
              {name: 'step4', run: 'ui-queries/final.cwl'}
            ]
          ]
        }
      };

      let result = completedUiStepsSelector(context);
      expect(result).toEqual([
        {step: {name: 'step3', run: 'ui-queries/other.cwl'}, index: 2}
      ]);
    });
  });

  describe('nextUiStepsSelector', () => {
    const context = {
      status: [
        [
          {
            name: 'step1',
            status: 'complete',
            downloads: {
              output: 'https://foo'
            },
            data: {
              output: [1, 2, 4]
            }
          },
          {
            name: 'step2',
            status: 'pending'
          },
          {
            name: 'step3',
            status: 'pending'
          },
          {
            name: 'step4',
            status: 'error'
          }
        ]
      ],
      pathIndex: 0,
      workflow: {
        steps: [
          [
            {name: 'step1', run: 'scripts/real.cwl', out: ['output']},
            {
              name: 'step2',
              run: 'ui-queries/something.cwl',
              in: [{source: ['step1', 'output']}]
            },
            {
              name: 'step3',
              run: 'ui-queries/other.cwl',
              in: [{source: ['step1', 'output']}]
            },
            {name: 'step4', run: 'ui-queries/final.cwl'}
          ]
        ]
      }
    };

    let result = nextUiStepsSelector(context);
    expect(result).toEqual([
      {
        name: 'step2',
        index: 1,
        run: 'ui-queries/something.cwl',
        in: [{source: ['step1', 'output']}]
      },
      {
        name: 'step3',
        index: 2,
        run: 'ui-queries/other.cwl',
        in: [{source: ['step1', 'output']}]
      }
    ]);
  });

  describe('errorStepsSelector', () => {
    const context = {
      status: [
        [
          {
            name: 'step1',
            status: 'pending'
          },
          {
            name: 'step2',
            status: 'complete'
          },
          {
            name: 'step3',
            status: 'complete'
          },
          {
            name: 'step4',
            status: 'error'
          }
        ]
      ],
      pathIndex: 0,
      workflow: {
        steps: [
          [
            {name: 'step1', run: 'ui-queries/something.cwl'},
            {name: 'step2', run: 'scripts/real.cwl'},
            {name: 'step3', run: 'ui-queries/other.cwl'},
            {name: 'step4', run: 'scripts/final.cwl'}
          ]
        ]
      }
    };

    let result = errorStepsSelector(context);
    expect(result).toEqual([
      {step: {name: 'step4', run: 'scripts/final.cwl'}, index: 3}
    ]);
  });
});
