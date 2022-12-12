import * as _ from 'lodash';

import {QuerySubclause} from '../../contexts/query/query_types';

export default class FilterOperator {
  isColumnFilter: boolean;
  subclause: QuerySubclause;

  static queryOperatorsByType: {[key: string]: {[key: string]: string}} = {
    base: {
      'Is present': '::has',
      'Is missing': '::lacks'
    },
    boolean: {
      'Is true': '::true',
      'Is false': '::false',
      'Is untrue': '::untrue'
    },
    number: {
      In: '::in',
      Equals: '::=',
      'Greater than': '::>',
      'Greater than or equals': '::>=',
      'Less than': '::<',
      'Less than or equals': '::<=',
      'Not equals': '::!=',
      'Not in': '::notin'
    },
    date: {
      Equals: '::=',
      'Greater than': '::>',
      'Greater than or equals': '::>=',
      'Less than': '::<',
      'Less than or equals': '::<=',
      'Not equals': '::!='
    },
    text: {
      In: '::in',
      Equals: '::equals',
      Contains: '::matches',
      Not: '::not',
      'Not in': '::notin',
      'Greater than': '::>',
      'Greater than or equals': '::>=',
      'Less than': '::<',
      'Less than or equals': '::<='
    }
  };

  static columnOptionsByType: {[key: string]: {[key: string]: string}} = {
    matrix: {
      Slice: '::slice'
    }
  };

  static terminalOperators: string[] = ['::true', '::false', '::untrue'];

  static terminalInvertOperators: string[] = ['::has', '::lacks'];

  static commaSeparatedOperators: string[] = ['::in', '::slice', '::notin'];

  static numericTypes: string[] = ['number', 'integer', 'float'];

  static allOperators(): string[] {
    return [
      ...new Set(
        Object.values(FilterOperator.queryOperatorsByType).reduce(
          (operators: string[], textOperatorMap) => {
            operators = operators.concat(Object.values(textOperatorMap));

            return operators;
          },
          []
        )
      )
    ];
  }

  constructor({
    subclause,
    isColumnFilter
  }: {
    subclause: QuerySubclause;
    isColumnFilter: boolean;
  }) {
    this.subclause = subclause;
    this.isColumnFilter = isColumnFilter;
  }

  hasOperand(): boolean {
    return !(
      FilterOperator.terminalOperators.includes(this.subclause.operator) ||
      FilterOperator.terminalInvertOperators.includes(this.subclause.operator)
    );
  }

  hasPrepopulatedOperandOptions(): boolean {
    return (
      'string' === this.subclause.attributeType &&
      '' !== this.subclause.attributeName
    );
  }

  attributeInputType(): string {
    switch (this.subclause.attributeType) {
      case 'string':
        return 'text';
      case 'date_time':
        return 'date';
      case 'integer':
      case 'float':
      case 'number':
        return 'number';
      case 'boolean':
        return 'boolean';
      case 'matrix':
        return 'matrix';
      case 'child':
      case 'table':
      case 'collection':
        return 'collection';
      default:
        return 'text';
    }
  }

  optionsForAttribute(): {[key: string]: string} {
    return this.isColumnFilter &&
      this.attributeInputType() in FilterOperator.columnOptionsByType
      ? FilterOperator.columnOptionsByType[this.attributeInputType()]
      : this.attrOptionsWithBaseOptions();
  }

  attrOptionsWithBaseOptions(): {[key: string]: string} {
    return {
      ...FilterOperator.queryOperatorsByType.base,
      ...(FilterOperator.queryOperatorsByType[this.attributeInputType()] || {})
    };
  }

  prettify(): string {
    return _.invert(this.optionsForAttribute())[this.subclause.operator];
  }

  magmify(newOperator: string): string {
    return this.optionsForAttribute()[newOperator];
  }

  options(): {[key: string]: string} {
    return this.optionsForAttribute();
  }

  formatOperand(operand: string): number | string {
    return operand;
  }
}
