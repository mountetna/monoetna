import {validator} from '../config-script';
import { Editor } from 'codemirror';

describe('Json Schema Validator', () => {
  let schema:any;
  let editor:any = {

  };
  beforeEach(() => {
    schema = {
      type: "object",
      properties: {
        foo: {type: "integer"},
        bar: { enum: [ "baz", "qux" ] }
      },
      required: ["foo"],
      additionalProperties: false
    }
  });

  it('ignores blank inputs', () => {
    const validate = validator(schema, editor as Editor)

    const text = '';
    const result = validate(text);

    expect(result).toEqual([]);
  });

  it('validates a correct json', () => {
    const validate = validator(schema, editor as Editor)

    const text = '{ "foo": 2 }';
    const result = validate(text);

    expect(result).toEqual([]);
  });

  it('returns errors for an incorrect json', () => {
    const validate = validator(schema, editor as Editor);

    const text = '{ "foo": "bar" }';
    const result = validate(text);

    expect(result[0].message).toEqual("/foo must be integer");
  });

  it('returns additionalProperty errors', () => {
    const validate = validator(schema, editor as Editor);

    const text = '{ "foo": 1, "baz": "qux" }';
    const result = validate(text);

    expect(result[0].message).toEqual("must NOT have additional properties baz");
  });

  it('returns enum errors', () => {
    const validate = validator(schema, editor as Editor);

    const text = '{ "foo": 1, "bar": "fred" }';
    const result = validate(text);

    expect(result[0].message).toEqual("/bar must be equal to one of the allowed values baz, qux");
  });
});
