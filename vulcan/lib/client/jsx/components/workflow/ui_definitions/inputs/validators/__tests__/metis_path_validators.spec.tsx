import {DataEnvelope, ValidationInputSpecification} from '../../input_types';
import {some} from '../../../../../../selectors/maybe';
import { MetisPathValidator, MetisFileValidator, MetisFolderValidator } from '../metis_path_validators';

describe('MetisPathValidator', () => {
  let input: ValidationInputSpecification<DataEnvelope<any>, DataEnvelope<any>>;

  beforeEach(() => {
    input = {
      value: [{
        bucket: '',
        path: '',
        type: null
      }],
      data: {}
    };
  });

  it('reports error for missing bucket', () => {
    expect(MetisPathValidator()(input).length > 0).toEqual(true);
  });

  it('allows path of \'\'', () => {
    input.value = some({bucket: 'data', path: '', type: 'folder'});
    expect(MetisPathValidator()(input)).toEqual([]);
  });

  it('allows target type = folder', () => {
    input.value = some({bucket: 'data', path: '', type: 'folder'});
    expect(MetisPathValidator()(input)).toEqual([]);
  });

  it('allows target type = file', () => {
    input.value = some({bucket: 'data', path: 'some_file.txt', type: 'file'});
    expect(MetisPathValidator()(input)).toEqual([]);
  });
});

describe('MetisFileValidator', () => {
  let input: ValidationInputSpecification<DataEnvelope<any>, DataEnvelope<any>>;

  beforeEach(() => {
    input = {
      value: [{
        bucket: '',
        path: '',
        type: 'file'
      }],
      data: {}
    };
  });

  it('reports error for missing bucket', () => {
    expect(MetisFileValidator()(input).length > 0).toEqual(true);
  });

  it('reports error for path of \'\'', () => {
    input.value = some({bucket: 'data', path: '', type: 'folder'});
    expect(MetisFileValidator()(input).length > 0).toEqual(true);
  });

  it('reports error for target type = folder', () => {
    input.value = some({bucket: 'data', path: 'subfolder', type: 'folder'});
    expect(MetisFileValidator()(input).length > 0).toEqual(true);
  });

  it('allows target type = file', () => {
    input.value = some({bucket: 'data', path: 'some_file.txt', type: 'file'});
    expect(MetisFileValidator()(input)).toEqual([]);
  });

  it('reports error with regex if regex given and not matched', () => {
    input.value = some({bucket: 'data', path: 'some_file.txt', type: 'file'});
    expect(MetisFileValidator('\\.[ct]sv$')(input).length > 0).toEqual(true);
    expect(MetisFileValidator('\\.[ct]sv$')(input)[0]).toMatch(/path must/);
    expect(MetisFileValidator('\\.[ct]sv$')(input)[0]).toMatch(/\\\.\[ct\]sv\$/);
  });

  it('reports error with descriptor and regex if regex given and not matched', () => {
    input.value = some({bucket: 'data', path: 'some_file.txt', type: 'file'});
    expect(MetisFileValidator('\\.[ct]sv$', 'csv or tsv file')(input).length > 0).toEqual(true);
    expect(MetisFileValidator('\\.[ct]sv$', 'csv or tsv file')(input)[0]).toMatch(/path should/);
    expect(MetisFileValidator('\\.[ct]sv$', 'csv or tsv file')(input)[0]).toMatch(/csv or tsv file/);
    expect(MetisFileValidator('\\.[ct]sv$')(input)[0]).toMatch(/\\\.\[ct\]sv\$/);
  });

  it('allows target of csv when given path_regex targeting csv or tsv', () => {
    input.value = some({bucket: 'data', path: 'some_file.csv', type: 'file'});
    expect(MetisFileValidator('\\.[ct]sv$', 'csv or tsv file')(input)).toEqual([]);
  });
});

describe('MetisFolderValidator', () => {
  let input: ValidationInputSpecification<DataEnvelope<any>, DataEnvelope<any>>;

  beforeEach(() => {
    input = {
      value: [{
        bucket: '',
        path: '',
        type: 'folder'
      }],
      data: {}
    };
  });

  it('reports error for missing bucket', () => {
    expect(MetisFolderValidator()(input).length > 0).toEqual(true);
  });

  it('allows path of \'\'', () => {
    input.value = some({bucket: 'data', path: '', type: 'folder'});
    expect(MetisFolderValidator()(input)).toEqual([]);
  });

  it('allows target type = folder', () => {
    input.value = some({bucket: 'data', path: '', type: 'folder'});
    expect(MetisFolderValidator()(input)).toEqual([]);
  });

  it('reports error for target type = file', () => {
    input.value = some({bucket: 'data', path: 'some_file.txt', type: 'file'});
    expect(MetisFolderValidator()(input).length > 0).toEqual(true);
  });
});
