import {createSelector} from 'reselect';

class SearchCache {
  constructor(search) {
    this.model_name = search.model_name;
    this.cached_pages = new Set(Object.keys(search.pages));
    this.current_page = search.current_page;
    this.record_names =
      search.current_page && search.pages[search.current_page]
        ? search.pages[search.current_page].record_names
        : null;
    this.page_size = search.page_size;
  }

  isCached(page) {
    return this.cached_pages.has(page);
  }
}

const selectSearchData = (state) => state.search;

export const selectSearchCache = createSelector(
  selectSearchData,
  (search) => new SearchCache(search)
);

export const selectSearchAttributeNames = createSelector(
  selectSearchData,
  (search) => search.attribute_names
);

export const selectSearchFilterString = createSelector(
  selectSearchData,
  (search) => search.filter_string
);

export const selectSearchFilterParams = createSelector(
  selectSearchData,
  (search) => search.filter_params
);

export const constructSingleFilterString = createSelector(
  selectSearchData,
  ({filter_params, filter_string}) => {
    // filter_params used for basic filter
    // filter_string used for advanced filter
    // Only one should be populated.

    if (filter_string) return filter_string;

    if (filter_params)
      return filter_params
        .map((param) => {
          return `${param.attribute}${param.operator}${param.value}`;
        })
        .join(' ');

    return '';
  }
);
