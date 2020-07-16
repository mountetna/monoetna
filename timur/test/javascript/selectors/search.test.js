import {
  selectSearchAttributeNames,
  constructSingleFilterString
} from '../../../lib/client/jsx/selectors/search';

describe('selectSearchAttributeNames', () => {
  it('returns the attribute_names from the search state', () => {
    let state = {
      search: {
        attribute_names: ['name', 'stats', 'species']
      }
    };

    let attribute_names = selectSearchAttributeNames(state);

    expect(attribute_names).toEqual(['name', 'stats', 'species']);
  });
});

describe('constructSingleFilterString', () => {
  it('returns the filter_string if set', () => {
    let state = {
      search: {
        filter_params: [{attribute: 'species', operator: '=', value: 'lion'}],
        filter_string: 'all'
      }
    };

    let filter_string = constructSingleFilterString(state);

    expect(filter_string).toEqual('all');
  });

  it('returns the filter_params if no filter_string', () => {
    let state = {
      search: {
        filter_params: [
          {attribute: 'species', operator: '=', value: 'lion'},
          {attribute: 'name', operator: '~', value: '/nem/'},
          {attribute: 'lives', operator: '<', value: '2'},
          {attribute: 'stats', operator: '>', value: '5'}
        ]
      }
    };

    let filter_string = constructSingleFilterString(state);

    expect(filter_string).toEqual('species=lion name~/nem/ lives<2 stats>5');
  });
});
