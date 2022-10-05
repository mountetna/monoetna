import {validator} from '../config-script';

import {Diagnostic} from '@codemirror/lint';

describe('Json Schema Validator', () => {
  let schema: any;
  let editor: any = {};

  class MockDoc {
    text: string;

    constructor(txt: string) {
      this.text = txt;
    }

    toString(): string {
      return this.text;
    }
  }

  class MockView {
    state: {doc: MockDoc};

    constructor(text: string) {
      this.state = {
        doc: new MockDoc(text)
      };
    }
  }

  beforeEach(() => {
    schema = {
      type: 'object',
      properties: {
        foo: {type: 'integer'},
        bar: {enum: ['baz', 'qux']}
      },
      required: ['foo'],
      additionalProperties: false
    };
  });

  it('ignores blank inputs', () => {
    const validate = validator(schema);

    editor = new MockView('');
    const result = validate(editor);

    expect(result).toEqual([]);
  });

  it('validates a correct json', () => {
    const validate = validator(schema);

    editor = new MockView('{ "foo": 2 }');
    const result = validate(editor);

    expect(result).toEqual([]);
  });

  it('returns errors for an incorrect json', () => {
    const validate = validator(schema);

    editor = new MockView('{ "foo": "bar" }');
    const result = validate(editor);

    expect((result as Diagnostic[])[0].message).toEqual('/foo must be integer');
  });

  it('returns additionalProperty errors', () => {
    const validate = validator(schema);

    editor = new MockView('{ "foo": 1, "baz": "qux" }');
    const result = validate(editor);

    expect((result as Diagnostic[])[0].message).toEqual(
      'must NOT have additional properties baz'
    );
  });

  it('returns enum errors', () => {
    const validate = validator(schema);

    editor = new MockView('{ "foo": 1, "bar": "fred" }');
    const result = validate(editor);

    expect((result as Diagnostic[])[0].message).toEqual(
      '/bar must be equal to one of the allowed values baz, qux'
    );
  });
});
