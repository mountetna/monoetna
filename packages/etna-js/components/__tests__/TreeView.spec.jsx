import React from 'react';
import {mount} from 'enzyme';
import TreeView, {getSelectedLeaves} from "../TreeView";

class TreeViewTest {
  constructor() {
    this.options = [];
    this.selectable = true;
    this.onSelectionsChange = jest.fn();
  }

  get mounted() {
    return (this._mounted = this._mounted || mount(<TreeView options={this.options} selectable={this.selectable}
                                                             onSelectionsChange={this.onSelectionsChange}/>));
  }


  findInput(path, node) {
    return this.mounted.find('input[type="checkbox"]').findWhere(
      i => i.props()['data-path'].join(',') === path.join(',') &&
        node === i.props()['data-node']);
  }

  get selected() {
    return this.mounted.find('input[type="checkbox"]').getElements().filter(e => e.props.checked).map(e => e.props['data-path'].concat([e.props['data-node']]))
  }

  get lastSelectionsByCallback() {
    const {calls} = this.onSelectionsChange.mock;
    return getSelectedLeaves(calls[calls.length - 1][0]);
  }
}

describe('TreeView', () => {
  let test;
  beforeEach(() => (test = new TreeViewTest()));

  describe('with some options', () => {
    beforeEach(() => {
      test.options = [
        ['Top Level 1', [
          ['Inner 1'],
          ['Inner 2', [
            ['Inner Inner 1'],
            ['Inner Inner 2'],
            ['Inner Inner 3'],
          ]],
          ['Inner 3'],
        ]],
        ['Top Level 2'],
      ];
    });

    it('defaults those options enabled', () => {
      expect(test.selected).toEqual([
        ['Top Level 1'],
        ['Top Level 1', 'Inner 1'],
        ['Top Level 1', 'Inner 2'],
        ['Top Level 1', 'Inner 2', 'Inner Inner 1'],
        ['Top Level 1', 'Inner 2', 'Inner Inner 2'],
        ['Top Level 1', 'Inner 2', 'Inner Inner 3'],
        ['Top Level 1', 'Inner 3'],
        ['Top Level 2'],
      ]);
    });

    describe('selecting a top level', () => {
      beforeEach(() => {
        test.findInput([], 'Top Level 1').simulate('change');
      })

      it('disables recursively child nodes', () => {
        expect(test.selected).toEqual([[
          'Top Level 2',
        ]]);

        expect(test.lastSelectionsByCallback).toEqual([
          'Top Level 2',
        ])
      });

      describe('selecting each leaf', () => {
        it('continually enables those leaves and eventually the parents again, too', () => {
          test.findInput(['Top Level 1', 'Inner 2'], 'Inner Inner 1').simulate('change');
          expect(test.selected).toEqual([
            [ 'Top Level 1', 'Inner 2', 'Inner Inner 1' ],
            [ 'Top Level 2', ]
          ]);
          expect(test.lastSelectionsByCallback).toEqual([
            'Inner Inner 1',
            'Top Level 2',
          ])

          test.findInput(['Top Level 1', 'Inner 2'], 'Inner Inner 2').simulate('change');
          expect(test.selected).toEqual([
            [ 'Top Level 1', 'Inner 2', 'Inner Inner 1' ],
            [ 'Top Level 1', 'Inner 2', 'Inner Inner 2' ],
            [ 'Top Level 2', ]
          ]);
          expect(test.lastSelectionsByCallback).toEqual([
            'Inner Inner 1',
            'Inner Inner 2',
            'Top Level 2',
          ])

          test.findInput(['Top Level 1'], 'Inner 1').simulate('change');
          test.findInput(['Top Level 1'], 'Inner 3').simulate('change');
          expect(test.selected).toEqual([
            [ 'Top Level 1', 'Inner 1' ],
            [ 'Top Level 1', 'Inner 2', 'Inner Inner 1' ],
            [ 'Top Level 1', 'Inner 2', 'Inner Inner 2' ],
            [ 'Top Level 1', 'Inner 3' ],
            [ 'Top Level 2', ]
          ]);

          expect(test.lastSelectionsByCallback).toEqual([
            'Inner 1',
            'Inner Inner 1',
            'Inner Inner 2',
            'Inner 3',
            'Top Level 2',
          ])

          test.findInput(['Top Level 1', 'Inner 2'], 'Inner Inner 3').simulate('change');
          expect(test.selected).toEqual([
            [ 'Top Level 1' ],
            [ 'Top Level 1', 'Inner 1' ],
            [ 'Top Level 1', 'Inner 2' ],
            [ 'Top Level 1', 'Inner 2', 'Inner Inner 1' ],
            [ 'Top Level 1', 'Inner 2', 'Inner Inner 2' ],
            [ 'Top Level 1', 'Inner 2', 'Inner Inner 3' ],
            [ 'Top Level 1', 'Inner 3' ],
            [ 'Top Level 2', ]
          ]);
          expect(test.lastSelectionsByCallback).toEqual([
            'Inner 1',
            'Inner Inner 1',
            'Inner Inner 2',
            'Inner Inner 3',
            'Inner 3',
            'Top Level 2',
          ])
        });
      })

      describe('and then re-selecting that top level', () => {
        beforeEach(() => {
          test.findInput([], 'Top Level 1').simulate('change');
        });

        it('enables recursively child nodes', () => {
          expect(test.selected).toEqual([
            ['Top Level 1'],
            ['Top Level 1', 'Inner 1'],
            ['Top Level 1', 'Inner 2'],
            ['Top Level 1', 'Inner 2', 'Inner Inner 1'],
            ['Top Level 1', 'Inner 2', 'Inner Inner 2'],
            ['Top Level 1', 'Inner 2', 'Inner Inner 3'],
            ['Top Level 1', 'Inner 3'],
            ['Top Level 2'],
          ]);

          expect(test.lastSelectionsByCallback).toEqual([
            'Inner Inner 1',
            'Inner Inner 2',
            'Inner Inner 3',
            'Inner 3',
            'Top Level 2',
          ])
        });
      })
    });
  })
})
