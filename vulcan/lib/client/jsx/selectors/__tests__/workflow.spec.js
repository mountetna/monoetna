import {
  completedUiStepsSelector,
  nextUiStepIndexSelector,
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

  describe('nextUiStepIndexSelector', () => {
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

    let result = nextUiStepIndexSelector(context);
    expect(result).toEqual(0);
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
